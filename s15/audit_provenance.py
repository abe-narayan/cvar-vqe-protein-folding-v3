"""S15 AUDIT / Part E -- a machine-checkable provenance record.

Emits `s15/results/audit_provenance.json`: for every load-bearing artefact, its path,
size, sha256 prefix, mtime, the code that WRITES it, whether that writer is config-keyed,
whether a partial run could overwrite a complete one, and which of the five pinned
constants depends on it.

    python -m s15.audit_provenance          # build the record
    python -m s15.audit_provenance check    # re-verify hashes against the record
"""
from __future__ import annotations
import os, sys, json, glob, hashlib, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)
REC = os.path.join(OUT, "audit_provenance.json")

CONSTANTS = ["shipped", "pool_best", "top75_best", "synthesis_fit", "n_zero_recall"]

#: node -> (kind, path-or-glob, written_by, config_keyed, atomic_write, feeds)
NODES = [
    ("pdbs/*.pdb", "PRIMARY INPUT", "pdbs/*.pdb", "-- deposited coordinates --",
     False, True, CONSTANTS),
    ("pdbs_ext/*.pdb", "PRIMARY INPUT", "pdbs_ext/*.pdb", "-- deposited coordinates --",
     False, True, CONSTANTS),
    ("prots/*.pdb", "PRIMARY INPUT (fragment library source)", "prots/*.pdb",
     "-- deposited coordinates --", False, True, CONSTANTS),
    ("peptide_db.npz", "DERIVED CACHE", "peptide_db.npz",
     "peptide_db.build() / core.data; scan of DIRS", False, False, CONSTANTS),
    ("peptide_clusters.json", "PINNED ARTEFACT", "peptide_clusters.json",
     "peptide_db.clusters(); PINNED by BRIEF", False, False, CONSTANTS),
    ("peptide_folds.json", "PINNED ARTEFACT", "peptide_folds.json",
     "peptide_db.folds(); PINNED by BRIEF", False, False, CONSTANTS),
    ("results/benchmark_manifest.json", "PINNED ARTEFACT", "results/benchmark_manifest.json",
     "PINNED by BRIEF; benchmark, its BYTES ARE READ here, to hash it; no benchmark content is loaded, parsed or reported — the output carries only a digest, a size and mtimes", False, False, []),
    ("s8/generate_univ/*.npz", "DERIVED CACHE", "s8/generate_univ/*.npz",
     "s8.generate.stage_univ", False, True,
     ["pool_best", "top75_best", "n_zero_recall", "shipped", "synthesis_fit"]),
    ("s8/generate_rama.npz", "DERIVED CACHE", "s8/generate_rama.npz",
     "s8.generate.stage_univ", False, False, []),
    ("s8/project_prior.json", "TRAINED DATA", "s8/project_prior.json",
     "s8.project.stage_prior (fold-disciplined counts)", False, False, ["synthesis_fit"]),
    ("distogram_models", "TRAINED MODEL", "distogram_models/*",
     "training pipeline (out of scope here)", False, False, ["shipped"]),
    ("s12/cache/disto_*.npz", "DERIVED CACHE", "s12/cache/disto_*.npz",
     "s12.instrument.distogram", False, False, ["shipped"]),
    ("bench_results/cache/<cfg_key>/*.json", "DERIVED CACHE, CONFIG-KEYED",
     "bench_results/cache/1fc9f2dcf489e2fb/*.json",
     "core.pipeline._worker_run -> _atomic_write_json", True, True,
     ["top75_best", "synthesis_fit", "n_zero_recall"]),
    ("s13/results/qarch_enum_*.npz", "DERIVED CACHE", "s13/results/qarch_enum_*.npz",
     "s13.qarch_enum.enumerate_target", False, False, []),
    ("s14/cache/obj_enum_*.npz", "DERIVED CACHE", "s14/cache/obj_enum_*.npz",
     "s14.obj_enum", False, False, []),
    ("s14/cache/ener_norm.json", "FITTED CONSTANT", "s14/cache/ener_norm.json",
     "s14.ener_norm", False, False, []),
    ("s12/results/*.json", "RESULT", "s12/results/*.json",
     "s12.instrument.write -- NO CONFIG KEY (documented hazard)", False, False, []),
    ("s13/results/*.json", "RESULT", "s13/results/*.json",
     "s13.qarch_lib.write -- NO CONFIG KEY, NO completeness flag", False, False, []),
    ("s14/results/*.json", "RESULT", "s14/results/*.json",
     "various -- check each", False, False, []),
    ("esm_cache.npz / s12/esm_bank.py", "DERIVED CACHE", "esm_cache.npz",
     "ESM embedding cache; 1.5 GB, do not load (BRIEF s3)", False, False, ["shipped"]),
]


def _hash(path, cap=None):
    h = hashlib.sha256(); n = 0
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            h.update(b); n += len(b)
            if cap and n > cap:
                break
    return h.hexdigest()[:16], os.path.getsize(path)


def build():
    rec = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "root": ROOT, "nodes": []}
    for name, kind, pat, writer, keyed, atomic, feeds in NODES:
        files = sorted(glob.glob(os.path.join(ROOT, pat)))
        node = {"node": name, "kind": kind, "pattern": pat, "written_by": writer,
                "config_keyed": keyed, "atomic_write": atomic,
                "feeds_pinned_constants": feeds, "n_files": len(files)}
        if files:
            node["total_bytes"] = int(sum(os.path.getsize(f) for f in files))
            node["mtime_min"] = time.strftime("%Y-%m-%d", time.localtime(
                min(os.path.getmtime(f) for f in files)))
            node["mtime_max"] = time.strftime("%Y-%m-%d", time.localtime(
                max(os.path.getmtime(f) for f in files)))
            # hash up to 8 representative files
            sample = files if len(files) <= 8 else files[:4] + files[-4:]
            node["sample_sha256_16"] = {os.path.relpath(f, ROOT).replace("\\", "/"): _hash(f)[0]
                                        for f in sample}
        else:
            node["MISSING"] = True
        rec["nodes"].append(node)

    # ---- edges: the five constants and exactly what each reads -------------
    rec["constant_dependencies"] = {
        "pool_best": {
            "value": 1.7108244199364904,
            "computed_as": "mean over 126 targets of min(rr[order[:500]])",
            "reads": ["s8/generate_univ/<pdb>.npz:rr", "s8/generate_univ/<pdb>.npz:order"],
            "recomputed_by_instrument": False,
            "reproduced_from_primary": True,
            "reproduction_evidence": "s15/results/audit_repro_univ.json -- 126/126 "
                                     "universes rebuilt bit-identically from pdbs/ + "
                                     "pdbs_ext/ + prots/, mean 1.7108244199364904",
            "oracle": True,
            "oracle_note": "rr is a NATIVE-derived label; pool_best is an ORACLE CEILING, "
                           "not a predictive number"},
        "top75_best": {
            "value": 2.3061526409453816,
            "computed_as": "mean over 126 of min(rr[pool][sub]) with sub from the "
                           "production cache",
            "reads": ["s8/generate_univ/<pdb>.npz:rr",
                      "bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json:sub"],
            "recomputed_by_instrument": False,
            "reproduced_from_primary": "PARTIAL -- 3/126 targets rerun through "
                                       "core.pipeline.run_target, sub bit-identical",
            "oracle": True},
        "shipped": {
            "value": 3.4540004952559396,
            "computed_as": "mean over 126 of rr[argmin(shipped_score(distogram, D))]",
            "reads": ["s12/cache/disto_<pdb>.npz:risk", "s8/generate_univ/<pdb>.npz:W",
                      "s8/generate_univ/<pdb>.npz:rr", "distogram_models (via ESM)"],
            "recomputed_by_instrument": "SCORE ONLY -- the distogram is a cache read",
            "reproduced_from_primary": True,
            "reproduction_evidence": "s15/results/audit_determinism_t{1,2,8}.json -- "
                                     "distogram recomputed for all 126 at three torch "
                                     "thread counts; 0 argmin flips; value identical",
            "oracle": False,
            "caveat": "the distogram tensor is NOT bit-stable across torch thread counts "
                      "(prob up to 2.0e-6, risk up to 1.3e-4); the constant is."},
        "synthesis_fit": {
            "value": 3.2040761603809194,
            "computed_as": "mean over 126 of CA-RMSD(fit_ca, nat_ca)",
            "reads": ["bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json:fit_ca",
                      "s8/generate_univ/<pdb>.npz:nat_ca"],
            "recomputed_by_instrument": "KABSCH ONLY -- fit_ca is a cache read",
            "reproduced_from_primary": "PARTIAL -- 3/126 targets rerun end to end through "
                                       "core.pipeline.run_target at PROD config; fit_ca "
                                       "BIT-IDENTICAL, cfg_key matches Config().key()",
            "oracle": False},
        "n_zero_recall": {
            "value": 18,
            "computed_as": "count of targets where no pool member within pool_best+BAND "
                           "(BAND=1.5 A) survives into sub",
            "reads": ["s8/generate_univ/<pdb>.npz:rr",
                      "bench_results/cache/1fc9f2dcf489e2fb/<pdb>.json:sub"],
            "recomputed_by_instrument": True,
            "reproduced_from_primary": "PARTIAL (same as top75_best)",
            "oracle": True,
            "caveat": "BAND = 1.5 A is an ARBITRARY constant; the count 18 is a function "
                      "of it (see s12/adv_fail18.py BANDS sweep)."},
    }
    with open(REC, "w") as fh:
        json.dump(rec, fh, indent=1)
    print("wrote", REC)
    for n in rec["nodes"]:
        print(f"  {n['node']:45s} {n['n_files']:6d} files "
              f"{'CONFIG-KEYED' if n['config_keyed'] else 'no key':13s} "
              f"{'atomic' if n['atomic_write'] else 'NON-ATOMIC'}")
    return rec


def check():
    with open(REC) as fh:
        rec = json.load(fh)
    bad = []
    for n in rec["nodes"]:
        for rel, want in (n.get("sample_sha256_16") or {}).items():
            p = os.path.join(ROOT, rel)
            if not os.path.exists(p):
                bad.append((rel, "MISSING")); continue
            got = _hash(p)[0]
            if got != want:
                bad.append((rel, f"{want} -> {got}"))
    print(json.dumps({"n_checked": sum(len(n.get("sample_sha256_16") or {}) for n in rec["nodes"]),
                      "n_changed": len(bad), "changed": bad[:40]}, indent=1))
    return bad


if __name__ == "__main__":
    (check if (len(sys.argv) > 1 and sys.argv[1] == "check") else build)()
