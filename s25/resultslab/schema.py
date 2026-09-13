"""s25/resultslab/schema.py -- THE RESULTS DATA MODEL.

Writes, into `results/summary/`:

    results.json      one record per (target, configuration), plus the run's provenance
    results.csv       the same records, flat, one row each
    leaderboard.json  one aggregate per configuration
    leaderboard.csv   the same, flat

ABSENT MEANS ABSENT.  Every (target, configuration) cell in the grid produces a record.  A
cell with no prediction gets `status: "absent"`, `rmsd: null` and an `absent_reason`.  It is
never filled with a mean, a fallback, or the incumbent's value.  `RECORD_KEYS` is the full
key set and `ST.save_atomic` checks it on every row, so a partial run cannot pass itself off
as a complete one -- the failure mode `s12.instrument.write` documents and that
`save_atomic` exists to close.

STATISTICS COME FROM `s24.stats_lib`.  Nothing here re-derives an SE, a CI or an MDE.
Each configuration is compared to `baseline` PAIRED over the targets where BOTH are present,
and the block reports, per `stats_lib.compare`: effect, SE, MDE = 2.8016 x SE, effect/MDE,
the iid CI, the FOLD-CLUSTERED CI beside it, per-fold effects, folds-same-sign, W/L/T,
worst degradation, Type-M and the verdict string.  A configuration whose fold CI includes
zero is reported as NOT MEASURED, in those words, in the leaderboard the professor sees.

BASIS.  Every record carries `basis`.  Records on different bases are never aggregated into
one leaderboard row: `build()` takes exactly one basis per call and stamps it on the file.
"""
from __future__ import annotations

import csv
import datetime
import json
import math
import os
import sys
from typing import Dict, List, Optional, Sequence

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                     # noqa: E402
from s24 import stats_lib as ST                     # noqa: E402
from s25.resultslab import exportlib as EX          # noqa: E402
from s25.resultslab import providers as PV          # noqa: E402

#: The full key set every record must carry.  `save_atomic` verifies it row by row.
RECORD_KEYS = (
    "target_id", "pdb_id", "sequence", "length", "fold",
    "configuration", "status", "absent_reason",
    "selector", "hamiltonians", "distogram_used", "cvar_alpha",
    "qubits", "ansatz_layers", "candidate_count", "retained_count",
    # PRIMARY -- the production result.  `basis` is REQUIRED and is never null on a present
    # record; `_require_basis` raises rather than writing a number a reader cannot place.
    "rmsd", "basis", "prediction_path",
    # SECONDARY -- the non-physical intermediate, carried on the SAME row so the pair cannot
    # be separated by anything downstream, including a careless CSV column selection.
    "rmsd_secondary", "basis_secondary", "prediction_path_secondary",
    "basis_delta", "ca_bond_mean", "ca_bond_min", "ca_bond_mean_secondary",
    "ca_bond_min_secondary", "ca_bond_mean_native", "pool_best",
    "native_path",
    "seed", "git_commit", "module_hash", "timestamp", "is_synthetic", "provenance",
    "source",
)

LEADERBOARD_KEYS = (
    "configuration", "status", "n_present", "n_absent", "n_targets", "basis",
    "mean", "median", "sd", "q1", "q3", "best", "worst",
    "basis_secondary", "mean_secondary", "median_secondary", "basis_delta_mean",
    "pool_gate", "pool_best_mean", "pool_margin", "n_pool_violations",
    "difficulty_gate", "corr_with_pool_best", "sd_per_target",
    "frac_under_2", "frac_under_3",
    "paired_n", "paired_effect", "paired_median_effect", "se", "mde", "effect_over_mde",
    "ci95_iid_lo", "ci95_iid_hi", "ci95_fold_lo", "ci95_fold_hi",
    "folds_same_sign", "n_folds", "wins", "losses", "ties",
    "worst_degradation", "type_m", "type_m_flag", "verdict", "is_synthetic",
)

SUMMARY = EX.SUMMARY_DIR


#: ONE timestamp for a whole build, not one per record.  Two builds of the same inputs then
#: differ only in this field and in `provenance`, so `--check-determinism` can compare the
#: scientific content byte for byte instead of having to trust that it is stable.
_RUN_TS = [datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")]


def _now() -> str:
    return _RUN_TS[0]


def new_run_timestamp() -> str:
    _RUN_TS[0] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    return _RUN_TS[0]


def _mod_hash() -> str:
    return EX.module_hash(os.path.join(HERE, "exportlib.py"),
                          os.path.join(HERE, "providers.py"),
                          os.path.join(HERE, "schema.py"))


def _q(x, p):
    return float(np.percentile(np.asarray(x, float), p))


# ------------------------------------------------------------------ record construction
def _absent_record(t, configuration, basis, reason, git, mh):
    return {
        "target_id": t["target_id"], "pdb_id": t["pdb_id"], "sequence": t["sequence"],
        "length": int(t["length"]), "fold": int(t["fold"]),
        "configuration": configuration, "status": "absent", "absent_reason": str(reason),
        "selector": None, "hamiltonians": None, "distogram_used": None, "cvar_alpha": None,
        "qubits": None, "ansatz_layers": None, "candidate_count": None,
        "retained_count": None,
        "rmsd": None, "basis": basis, "prediction_path": None,
        "rmsd_secondary": None, "basis_secondary": EX.SECONDARY_BASIS,
        "prediction_path_secondary": None, "basis_delta": None,
        "ca_bond_mean": None, "ca_bond_min": None,
        "ca_bond_mean_secondary": None, "ca_bond_min_secondary": None,
        "ca_bond_mean_native": None, "pool_best": None,
        "native_path": None,
        "seed": None, "git_commit": git, "module_hash": mh, "timestamp": _now(),
        "is_synthetic": PV.UNKNOWN, "provenance": EX.PROV_UNKNOWN, "source": None,
    }


def _rel(p: Optional[str]) -> Optional[str]:
    return None if p is None else os.path.relpath(p, ROOT).replace("\\", "/")


def collect(configurations: Optional[Sequence[str]] = None,
            pdbs: Optional[Sequence[str]] = None,
            chain_lam: float = EX.PRODUCTION_LAM,
            verify: bool = True,
            require_genuine_header: bool = True,
            sandbox: bool = False,
            verbose: bool = True) -> List[Dict[str, object]]:
    """Walk the (target x configuration) grid, export BOTH bases, and build the records.

    There is no `basis` argument any more, and that is deliberate.  s25 L4/L8 rules that the
    built chain is the production result and the point cloud is a labelled non-physical
    intermediate reported ALWAYS PAIRED with it.  A `basis` switch would let a caller obtain
    one number without the other, which is exactly the failure the ruling exists to prevent.
    `EX.export_pair` emits both from one emitted point cloud and both land on the same row.

    `verify=True` re-reads every written file and asserts the round trip on BOTH bases, so
    the RMSDs in the record are the RMSDs of the bytes on disk.
    """
    cfgs = list(configurations) if configurations is not None else PV.registered()
    rows = all_rows = EX.target_rows()
    if pdbs is not None:
        want = {str(p).upper() for p in pdbs}
        rows = [r for r in all_rows if r["pdb_id"] in want]
    prov = ST.provenance(os.path.join(HERE, "schema.py"))
    git, mh = prov.get("git_commit"), _mod_hash()
    pb = EX.pool_best()

    # Natives first -- the viewer's overlay and the round-trip check both need them.  In a
    # SANDBOX build they go to the fixture directory too: results/ carries no data at all
    # until the frozen build, not even referents.
    for t in rows:
        EX.export_native(t["pdb_id"], sandbox=sandbox)

    out: List[Dict[str, object]] = []
    for cfg in cfgs:
        f = PV.get(cfg)
        smeta = PV.static_meta(cfg)
        n_ok = 0
        for t in rows:
            pdb = t["pdb_id"]
            try:
                got = f(pdb)
            except Exception as e:                                  # a provider failure is ABSENT
                out.append(_absent_record(t, cfg, EX.PRIMARY_BASIS,
                                          "provider raised: %s" % e, git, mh))
                continue
            if got is None or got.get("ca") is None:
                out.append(_absent_record(t, cfg, EX.PRIMARY_BASIS,
                                          "no prediction for this target", git, mh))
                continue
            meta = dict(smeta)
            meta.update({k: v for k, v in got.items() if k != "ca"})
            ca = np.asarray(got["ca"], float)

            prov = _provenance_of(meta)
            pair = EX.export_pair(pdb, cfg, ca, lam=chain_lam, provenance=prov)
            if verify:
                _verify_pair(pdb, cfg, pair, sandbox=sandbox)
                # GATE (c), read back off disk, not from the dict we meant to write
                allow = (EX.PROV_GENUINE,) if require_genuine_header else (
                    EX.PROV_GENUINE, EX.PROV_SYNTHETIC, EX.PROV_UNKNOWN)
                for which in ("primary", "secondary"):
                    EX.verify_header(pair[which]["path"], allow_provenance=allow)
                EX.verify_header(EX.native_path(pdb, sandbox=sandbox))

            rec = {
                "target_id": t["target_id"], "pdb_id": pdb, "sequence": t["sequence"],
                "length": int(t["length"]), "fold": int(t["fold"]),
                "configuration": cfg, "status": "present", "absent_reason": None,
                "selector": meta.get("selector"),
                "hamiltonians": meta.get("hamiltonians"),
                "distogram_used": meta.get("distogram_used"),
                "cvar_alpha": meta.get("cvar_alpha"),
                "qubits": meta.get("qubits"),
                "ansatz_layers": meta.get("ansatz_layers"),
                "candidate_count": meta.get("candidate_count"),
                "retained_count": meta.get("retained_count"),
                "rmsd": float(pair["rmsd"]), "basis": EX.PRIMARY_BASIS,
                "prediction_path": _rel(pair["primary"]["path"]),
                "rmsd_secondary": float(pair["rmsd_secondary"]),
                "basis_secondary": EX.SECONDARY_BASIS,
                "prediction_path_secondary": _rel(pair["secondary"]["path"]),
                "basis_delta": float(pair["delta"]),
                "ca_bond_mean": pair["primary"]["ca_bond"]["ca_bond_mean"],
                "ca_bond_min": pair["primary"]["ca_bond"]["ca_bond_min"],
                "ca_bond_mean_secondary": pair["secondary"]["ca_bond"]["ca_bond_mean"],
                "ca_bond_min_secondary": pair["secondary"]["ca_bond"]["ca_bond_min"],
                "ca_bond_mean_native": EX.ca_bond_stats(EX.native_ca(pdb))["ca_bond_mean"],
                "pool_best": pb.get(pdb),
                "native_path": _rel(EX.native_path(pdb, sandbox=sandbox)),
                "seed": meta.get("seed"), "git_commit": git, "module_hash": mh,
                "timestamp": _now(),
                "is_synthetic": meta.get("is_synthetic", PV.UNKNOWN),
                "provenance": prov,
                "source": meta.get("source"),
            }
            _require_basis(rec)
            out.append(rec)
            n_ok += 1
        if verbose:
            print("  %-24s %3d/%3d present" % (cfg, n_ok, len(rows)))
    return out


def _provenance_of(meta: Dict[str, object]) -> str:
    """Map a producer's declaration onto the header vocabulary.  Undeclared is UNKNOWN."""
    got = meta.get("is_synthetic", PV.UNKNOWN)
    if got is True:
        return EX.PROV_SYNTHETIC
    if got is False:
        return EX.PROV_GENUINE
    return EX.PROV_UNKNOWN


def _verify_pair(pdb, cfg, pair, sandbox: bool = False) -> None:
    """Round-trip both exported files against the exported native.  Raises, never warns."""
    nat_file = EX.read_pdb(EX.native_path(pdb, sandbox=sandbox))
    for which, key in (("primary", "rmsd"), ("secondary", "rmsd_secondary")):
        got = EX.read_pdb(pair[which]["path"])
        chk = I.ca_rmsd(got["CA"], nat_file["CA"])
        if abs(chk - float(pair[key])) >= EX.RMSD_TOL:
            raise AssertionError(
                "ROUND-TRIP FAILED for %s/%s [%s]: file-vs-file %.6f, exported header %.6f"
                % (pdb, cfg, pair[which]["basis"], chk, float(pair[key])))


def _require_basis(rec: Dict[str, object]) -> None:
    """A present record without a basis is not a record.  L4 exists because a number travelled
    without one; nothing here may emit an RMSD a reader has to guess the basis of."""
    if rec["status"] != "present":
        return
    for r_key, b_key in (("rmsd", "basis"), ("rmsd_secondary", "basis_secondary")):
        if rec.get(r_key) is None:
            raise ValueError("%s/%s: %s is null on a present record"
                             % (rec["target_id"], rec["configuration"], r_key))
        if not rec.get(b_key):
            raise ValueError("%s/%s: %s carries no %s -- refusing to emit an unplaceable RMSD"
                             % (rec["target_id"], rec["configuration"], r_key, b_key))


# ------------------------------------------------------------------ leaderboard
def leaderboard(records: Sequence[Dict[str, object]], baseline: str,
                n_targets: int, hard_gate: bool = True) -> List[Dict[str, object]]:
    """Aggregate per configuration: paired statistics from `s24.stats_lib`, and THE GATE.

    Every configuration is run through `EX.pool_oracle_gate` on the PRIMARY basis before any
    of its numbers are written.  With `hard_gate=True` (the frozen build) a configuration
    that scores below the ORACLE best member of its own K=500 candidate pool raises
    `PoolOracleViolation` and the build stops.  Nothing has to remember to set a flag.
    """
    basis = EX.PRIMARY_BASIS
    by_cfg: Dict[str, Dict[str, Dict[str, object]]] = {}
    for r in records:
        by_cfg.setdefault(r["configuration"], {})[r["pdb_id"]] = r

    base = by_cfg.get(baseline, {})
    rows: List[Dict[str, object]] = []
    order = [c for c in EX.CONFIGURATIONS if c in by_cfg] + [c for c in EX.FIXTURES if c in by_cfg]
    for cfg in order:
        d = by_cfg[cfg]
        present = {p: r for p, r in d.items() if r["status"] == "present"}
        vals = np.array([r["rmsd"] for r in present.values()], float)
        sec = np.array([r["rmsd_secondary"] for r in present.values()], float)
        row = {k: None for k in LEADERBOARD_KEYS}
        prov = {r.get("is_synthetic") for r in present.values()}
        row.update({
            "configuration": cfg, "basis": basis, "basis_secondary": EX.SECONDARY_BASIS,
            "n_present": int(len(present)), "n_absent": int(len(d) - len(present)),
            "n_targets": int(n_targets),
            "status": "complete" if len(present) == n_targets else "partial",
            "is_synthetic": (True if True in prov else
                             (PV.UNKNOWN if PV.UNKNOWN in prov else False)),
        })
        if len(vals):
            per = {p: r["rmsd"] for p, r in present.items()}
            gate = EX.pool_oracle_gate(per, cfg, hard=hard_gate)
            diff = EX.difficulty_gate(per, cfg, hard=hard_gate)
            row.update({
                "mean": float(vals.mean()), "median": float(np.median(vals)),
                "sd": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0,
                "q1": _q(vals, 25), "q3": _q(vals, 75),
                "best": float(vals.min()), "worst": float(vals.max()),
                "mean_secondary": float(sec.mean()), "median_secondary": float(np.median(sec)),
                "basis_delta_mean": float((vals - sec).mean()),
                "pool_gate": gate["status"], "pool_best_mean": gate["pool_best_mean"],
                "pool_margin": gate["margin"],
                "n_pool_violations": gate["n_per_target_violations"],
                "difficulty_gate": diff["status"],
                "corr_with_pool_best": diff.get("corr_with_pool_best"),
                "sd_per_target": diff.get("sd"),
                "frac_under_2": float((vals < 2.0).mean()),
                "frac_under_3": float((vals < 3.0).mean()),
            })
            row["_gate"] = gate
            row["_difficulty"] = diff
        # paired against the baseline, over targets present in BOTH
        shared = sorted(p for p in present if p in base and base[p]["status"] == "present")
        if cfg != baseline and len(shared) >= 8:
            a = np.array([present[p]["rmsd"] for p in shared], float)
            b = np.array([base[p]["rmsd"] for p in shared], float)
            folds = np.array([present[p]["fold"] for p in shared], int)
            c = ST.compare(a, b, folds, names=shared,
                           label="%s vs %s [%s]" % (cfg, baseline, basis))
            row.update({
                "paired_n": c["n"], "paired_effect": c["effect"],
                "paired_median_effect": c["median_effect"],
                "se": c["se"], "mde": c["mde"], "effect_over_mde": c["effect_over_mde"],
                "ci95_iid_lo": c["ci95_iid"][0], "ci95_iid_hi": c["ci95_iid"][1],
                "ci95_fold_lo": c["ci95_fold"][0] if c["ci95_fold"] else None,
                "ci95_fold_hi": c["ci95_fold"][1] if c["ci95_fold"] else None,
                "folds_same_sign": c.get("folds_same_sign"), "n_folds": c.get("n_folds"),
                "wins": c["n_better"], "losses": c["n_worse"], "ties": c["n_tied"],
                "worst_degradation": c["worst_degradation"],
                "type_m": c["type_m"], "type_m_flag": c["type_m_flag"],
                "verdict": c["verdict"],
            })
            row["_compare"] = c
        elif cfg == baseline:
            row["verdict"] = "BASELINE"
        else:
            row["verdict"] = "NOT MEASURED (fewer than 8 paired targets)"
        rows.append(row)
    rows.sort(key=lambda r: (r["mean"] if r["mean"] is not None else math.inf))
    return rows


# ------------------------------------------------------------------ writing
def _flat(v):
    if isinstance(v, (list, tuple)):
        return "|".join(str(x) for x in v)
    if isinstance(v, bool):
        return "true" if v else "false"
    return "" if v is None else v


def _write_csv(path: str, rows: Sequence[Dict[str, object]], keys: Sequence[str]) -> str:
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([_flat(r.get(k)) for k in keys])
    os.replace(tmp, path)
    return path


def build(baseline: str = "production",
          configurations: Optional[Sequence[str]] = None,
          pdbs: Optional[Sequence[str]] = None, label: str = "",
          chain_lam: float = EX.PRODUCTION_LAM, out_dir: Optional[str] = None,
          require_provenance: bool = True, hard_gate: bool = True,
          allow_synthetic: bool = False,
          verbose: bool = True) -> Dict[str, object]:
    """Run the grid, write the four summary artefacts, return the payload.

    THREE REFUSALS, in the order they fire:

    * `require_provenance` -- a present record whose `is_synthetic` is UNKNOWN (a producer
      that did not declare it) stops the build.  Absent provenance is not genuine provenance.
    * `allow_synthetic` -- a build that is not explicitly a selftest refuses to write a
      record declared synthetic.  Together with `providers.register` refusing to bind a
      synthetic provider to a real configuration name, a synthetic number cannot reach a
      shipped artefact by any route.
    * `hard_gate` -- THE POOL-ORACLE GATE, applied in `leaderboard`.  A configuration
      scoring below the ORACLE best member of its own K=500 candidate pool raises and the
      build stops.  It is native-free and machine-checkable, it needs nobody to remember a
      flag, and it subsumes the other two: it fires on leakage nobody anticipated.
    """
    new_run_timestamp()
    cfgs = list(configurations) if configurations is not None else PV.registered()
    if baseline not in cfgs:
        raise KeyError("baseline %r is not among the configurations being built" % baseline)
    summary_dir = out_dir or SUMMARY
    os.makedirs(summary_dir, exist_ok=True)
    tr = EX.target_rows()
    n_targets = len(tr) if pdbs is None else len(pdbs)

    records = collect(cfgs, pdbs=pdbs, chain_lam=chain_lam, verbose=verbose,
                      require_genuine_header=not allow_synthetic, sandbox=allow_synthetic)

    n_present = sum(1 for r in records if r["status"] == "present")
    unknown = [r for r in records
               if r["status"] == "present" and r.get("is_synthetic") == PV.UNKNOWN]
    if require_provenance and unknown:
        raise ValueError(
            "PROVENANCE UNDECLARED on %d of %d present records (e.g. %s/%s from %r). A "
            "producer that does not say whether its structures are real is UNKNOWN, never "
            "genuine. Declare `is_synthetic` in the spec meta block or in the emitted file."
            % (len(unknown), n_present, unknown[0]["target_id"],
               unknown[0]["configuration"], unknown[0]["source"]))
    synth = [r for r in records if r.get("is_synthetic") is True]
    if synth and not allow_synthetic:
        raise ValueError(
            "REFUSED: %d records are declared SYNTHETIC and this is not a selftest build. "
            "Synthetic structures may not reach results/. Use the selftest mode, which "
            "writes to its own sandbox under test-fixture labels." % len(synth))

    lb = leaderboard(records, baseline=baseline, n_targets=n_targets, hard_gate=hard_gate)
    synthetic = bool(synth)

    payload = {
        "label": label or ("PROVISIONAL -- synthetic rows present" if synthetic else ""),
        "basis": EX.PRIMARY_BASIS,
        "basis_note": ("BUILT CHAIN (lam=" + ("%.2f" % chain_lam) + "): the emitted point "
                       "cloud projected to the nearest ideal-geometry backbone. THE "
                       "PRODUCTION RESULT -- a real molecule that can be written as a valid "
                       "PDB and round-tripped as one. Always reported paired with the point "
                       "cloud, which reads about 0.166 A better and is a NON-PHYSICAL "
                       "INTERMEDIATE: mean virtual Ca-Ca bond 2.9614 A against a native "
                       "3.8122 A, 22.3 per cent contracted, worst single bond 0.649 A."),
        "basis_secondary": EX.SECONDARY_BASIS,
        "basis_secondary_note": ("POINT CLOUD: the emitted coordinate average. Retained and "
                                 "reported everywhere, never the headline. Not a molecule."),
        "chain_lam": float(chain_lam),
        "baseline": baseline,
        "configurations": cfgs,
        "n_targets": n_targets,
        "target_map_digest": EX.target_map()["digest"],
        "contains_synthetic": synthetic,
        "status": "PROVISIONAL" if synthetic else "GENERATED",
        "mde_rule": "MDE = 2.8016 * SE, per comparison (s24/stats_lib)",
        "pool_gate_rule": ("no configuration may score below the ORACLE best member of its "
                           "own K=500 candidate pool (mean " + ("%.4f" % EX.POOL_BEST_MEAN)
                           + " A over the 126 dev targets). PASS = no target below its own "
                           "pool best; WARN = one or more individual targets below their own "
                           "pool best, which is not a defect: the emitted chain is a "
                           "coordinate average projected onto ideal geometry, not a pool "
                           "member, so on a single target it can land nearer the native than "
                           "any one window; FAIL = the mean is below the pool-best mean, "
                           "which no real method can do"),
        "pool_gate": {r["configuration"]: r.get("_gate") for r in lb},
        "difficulty_gate": {r["configuration"]: r.get("_difficulty") for r in lb},
        "difficulty_gate_rule": (
            "corr(per-target RMSD, that target's ORACLE pool_best) >= %.2f AND "
            "sd(per-target RMSD) >= %.2f; incumbent +0.6988/1.6466 (point cloud), "
            "+0.7157/1.7361 (production chain); synthetic arms -0.12..+0.08 / 0.20..0.47"
            % (EX.DIFFICULTY_CORR_FLOOR, EX.DIFFICULTY_SD_FLOOR)),
        "records": records,
        "leaderboard": [{k: v for k, v in r.items() if not k.startswith("_")} for r in lb],
    }
    rp = os.path.join(summary_dir, "results.json")
    ST.save_atomic(rp, payload, complete_keys=RECORD_KEYS, rows=records,
                   n_expected=n_targets * len(cfgs), module_file=os.path.join(HERE, "schema.py"))
    _write_csv(os.path.join(summary_dir, "results.csv"), records, RECORD_KEYS)

    lp = os.path.join(summary_dir, "leaderboard.json")
    ST.save_atomic(lp, {"basis": EX.PRIMARY_BASIS, "basis_secondary": EX.SECONDARY_BASIS,
                        "baseline": baseline, "n_targets": n_targets,
                        "status": payload["status"], "rows": payload["leaderboard"]},
                   complete_keys=LEADERBOARD_KEYS, rows=payload["leaderboard"],
                   n_expected=len(cfgs), module_file=os.path.join(HERE, "schema.py"))
    _write_csv(os.path.join(summary_dir, "leaderboard.csv"), payload["leaderboard"],
               LEADERBOARD_KEYS)

    if verbose:
        print(fmt_leaderboard(payload))
    return payload


def fmt_leaderboard(payload: Dict[str, object]) -> str:
    L = ["", "  LEADERBOARD   baseline=%s   status=%s"
         % (payload["baseline"], payload["status"]),
         "  PRIMARY %s (production, a real chain)  |  secondary %s (NON-PHYSICAL)"
         % (payload["basis"], payload["basis_secondary"]),
         "  %-22s %5s %8s %8s %9s %9s %9s %5s %5s %7s  %s"
         % ("configuration", "n", "mean", "median", "pt-cloud", "effect", "MDE", "pool",
            "diff", "corr", "verdict")]
    for r in payload["leaderboard"]:
        L.append("  %-22s %5s %8s %8s %9s %9s %9s %5s %5s %7s  %s"
                 % (r["configuration"],
                    "%d/%d" % (r["n_present"], r["n_targets"]),
                    "--" if r["mean"] is None else "%8.4f" % r["mean"],
                    "--" if r["median"] is None else "%8.4f" % r["median"],
                    "--" if r["mean_secondary"] is None else "%9.4f" % r["mean_secondary"],
                    "--" if r["paired_effect"] is None else "%+9.4f" % r["paired_effect"],
                    "--" if r["mde"] is None else "%9.4f" % r["mde"],
                    (r["pool_gate"] or "--")[:5],
                    (r["difficulty_gate"] or "--")[:5],
                    "--" if r["corr_with_pool_best"] is None else "%+7.3f" % r["corr_with_pool_best"],
                    r["verdict"] or ""))
    warn = [r for r in payload["leaderboard"] if r.get("pool_gate") == "WARN"]
    if warn:
        L.append("  pool=WARN is per-target and is not a defect: %s beat the ORACLE best member "
                 "of the target's own pool on that many targets. The emitted chain is a "
                 "coordinate average projected onto ideal geometry, not a pool member, so on a "
                 "single target it can land nearer the native than any one window. The release "
                 "condition is the aggregate (mean >= pool-best mean), which every row passes."
                 % ", ".join("%s (%d target%s)" % (r["configuration"], r["n_pool_violations"] or 0,
                                                    "" if (r["n_pool_violations"] or 0) == 1 else "s")
                             for r in warn))
    if payload["contains_synthetic"]:
        L.append("  *** PROVISIONAL: synthetic rows present.  Not a result. ***")
    return "\n".join(L)
