"""s25/resultslab/build.py -- THE ONE COMMAND.  Registers providers, runs the grid, renders.

    python -m s25.resultslab.build --mode selftest --limit 12   # sandbox smoke test
    python -m s25.resultslab.build --mode selftest              # sandbox, all 126
    python -m s25.resultslab.build --mode frozen --spec s25/results/frozen_configs.json

PHASE ORDER (BRIEF SS7, and the coordinator's ruling of this session).  The architecture is
not frozen.  Therefore:

  * `results/` SHIPS WITH NO DATA until the freeze signal.  `--mode frozen` is the only mode
    that writes there, and it requires a spec file that does not yet exist.
  * NO structure may be exported under a real configuration name before the signal.
  * `--mode selftest` proves the export -> schema -> leaderboard -> render path end to end
    WITHOUT touching any of that: it writes to its own sandbox at
    `s25/resultslab/_selftest/`, and it registers its arms under TEST FIXTURE labels
    (`_incumbent`, `_zero`, `_noise`) which `exportlib` quarantines in
    the SYSTEM TEMP DIRECTORY, and which `providers.register` is the only path to.

The earlier `--mode proof` registered native-plus-noise under six of the seven REAL
configuration names and rendered them into `results/site`.  Nothing was wrong with the
numbers -- they were flagged synthetic everywhere -- but the only thing standing between a
synthetic structure and a page reading "legacy: 2.10 A" was a flag a reader had to notice.
That mode is gone.  `providers.register` now refuses the binding outright, and
`schema.build` refuses to write a synthetic record unless `allow_synthetic` is set, which
only this module's selftest path sets.

THE SPEC FILE, for `--mode frozen`.  It, and not this module, decides what a configuration
IS.  Nothing in this package changes when the architecture freezes:

    {
      "baseline": "production",
      "label": "S25 frozen architecture",
      "architecture": {"selector": "...", "...": "..."},
      "configurations": {
        "legacy": {"path": "s25/results/emit_legacy.json",
                   "meta": {"selector": "cvar_vqe", "hamiltonians": ["legacy"],
                            "distogram_used": false, "cvar_alpha": 0.15,
                            "qubits": 9, "ansatz_layers": 3, "seed": 0,
                            "is_synthetic": false}},
        "production": {"path": "..."}
      }
    }

Each `path` is a JSON of `{"<PDB id>": [[x,y,z], ...]}` -- the EMITTED POINT CLOUD, in the
same residue order as the target's sequence.  The lab projects it to the built chain itself,
through the real stage-3b path, so every configuration is chained by the same operator.
A target missing from the file is ABSENT and is written as absent.
`meta.is_synthetic` is REQUIRED: absent provenance is UNKNOWN and fails the build.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s25.resultslab import exportlib as EX          # noqa: E402
from s25.resultslab import providers as PV          # noqa: E402
from s25.resultslab import schema as SC             # noqa: E402
from s25.resultslab import site as SITE             # noqa: E402

SELFTEST_DIR = os.path.join(HERE, "_selftest")


def register_selftest() -> None:
    """Three FIXTURE arms.  None of them may carry a real configuration name; `register`
    enforces that for the two synthetic ones, and `_incumbent` is a replay, not a result."""
    PV.clear()
    PV.register("_incumbent", PV.incumbent(),
                selector="distogram_bayes_risk_argmin (classical control)",
                hamiltonians=[], distogram_used=True)
    PV.register("_zero", PV.synthetic(sigma=0.9, seed=1, tag="zero"))
    PV.register("_noise", PV.synthetic(sigma=2.4, seed=2, tag="noise"))


def register_frozen(spec_path: str) -> dict:
    with open(spec_path) as fh:
        spec = json.load(fh)
    PV.clear()
    for cfg, d in spec["configurations"].items():
        p = d["path"]
        p = p if os.path.isabs(p) else os.path.join(ROOT, p)
        loader = PV.from_npz if p.endswith(".npz") else PV.from_json
        meta = dict(d.get("meta", {}))
        #: the meta must reach the LOADER as well as the registry.  A provider's per-call
        #: `setdefault("is_synthetic", UNKNOWN)` is applied AFTER the static meta is merged in
        #: `schema.collect`, so a declaration passed only to `register` is silently overridden
        #: by the loader's own UNKNOWN default and the provenance gate then refuses genuine data.
        PV.register(cfg, loader(p, **meta), **meta)
    return spec


def _assert_gates_discriminate(payload) -> None:
    """The selftest's real assertion: the gates must SEPARATE the replay from the noise.

    A gate that never fires is not evidence.  This is what makes the selftest a test of the
    safeguards rather than a demonstration of the plumbing.
    """
    lb = {r["configuration"]: r for r in payload["leaderboard"]}
    print("")
    print("  GATE DEMONSTRATION (selftest, gates run SOFT so the table can exist)")
    print("  %-12s %-9s %-6s %-6s %9s %8s  %s"
          % ("label", "provenance", "pool", "diff", "corr", "sd", "expected"))
    for cfg, want in (("_incumbent", "PASS"), ("_zero", "FAIL"), ("_noise", "FAIL")):
        r = lb.get(cfg)
        if r is None:
            continue
        got = "PASS" if r["difficulty_gate"] == "PASS" else "FAIL"
        print("  %-12s %-9s %-6s %-6s %+9.4f %8.4f  %s %s"
              % (cfg, str(r["is_synthetic"]), r["pool_gate"], r["difficulty_gate"],
                 r["corr_with_pool_best"], r["sd_per_target"], want,
                 "OK" if got == want else "*** MISMATCH ***"))
        assert got == want, (
            "SELFTEST FAILED: the difficulty gate said %s for %r, expected %s. The gates are "
            "the release condition; if they stop discriminating here they are not protecting "
            "the frozen build either." % (got, cfg, want))
    print("  gates discriminate as designed")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build the Sprint-25 results laboratory.")
    ap.add_argument("--mode", choices=("selftest", "frozen"), default="selftest")
    ap.add_argument("--spec", default=None, help="config spec JSON (required for --mode frozen)")
    ap.add_argument("--chain-lam", type=float, default=EX.PRODUCTION_LAM,
                    help="torsion-prior weight for the built chain (0.3 = SYNTHESIS stage 3, "
                         "the shipped arm at 3.2148 A; 0.0 is the fit arm at 3.2041 A)")
    ap.add_argument("--baseline", default=None)
    ap.add_argument("--limit", type=int, default=None, help="first N targets only (smoke test)")
    ap.add_argument("--no-site", action="store_true")
    args = ap.parse_args(argv)

    EX.build_target_map()
    arch = None
    if args.mode == "selftest":
        register_selftest()
        baseline = args.baseline or "_incumbent"
        label = ("SELFTEST BUILD -- fixture labels only, sandboxed. Contains SYNTHETIC "
                 "structures. NOT A RESULT and not in results/.")
        out_dir = os.path.join(SELFTEST_DIR, "summary")
        site_dir = os.path.join(SELFTEST_DIR, "site")
        allow_synth = True
        # SOFT gates in the selftest -- on purpose.  The point of this mode is to SHOW the
        # gates firing on the synthetic arms; a hard gate would abort the run before the
        # table that demonstrates it exists.  The frozen build runs them hard.
        hard = False
    else:
        if not args.spec:
            ap.error("--mode frozen requires --spec (and the architecture must be frozen)")
        spec = register_frozen(args.spec)
        baseline = args.baseline or spec.get("baseline", "production")
        arch = spec.get("architecture")
        label = spec.get("label", "")
        out_dir, site_dir = EX.SUMMARY_DIR, EX.SITE_DIR
        allow_synth = False
        hard = True

    pdbs = None
    if args.limit:
        pdbs = [t["pdb_id"] for t in EX.target_rows()[:args.limit]]

    print("  mode=%s  labels: %s" % (args.mode, ", ".join(PV.registered())))
    payload = SC.build(baseline=baseline, pdbs=pdbs, label=label, chain_lam=args.chain_lam,
                       out_dir=out_dir, allow_synthetic=allow_synth,
                       require_provenance=True, hard_gate=hard)
    if args.mode == "selftest":
        _assert_gates_discriminate(payload)
    if not args.no_site:
        SITE.render(payload, out_dir=site_dir, architecture=arch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
