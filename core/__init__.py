"""`core` -- the consolidated production path of the peptide structure predictor.

Consolidated front half of the pipeline: data, geometry, prediction, cache; consolidated
back half: energy, AMBER, quantum; and `core.pipeline`, which is the whole of it.  Each
operation has exactly one implementation.  Importing a submodule is cheap -- nothing loads
a database, a model or the ESM bank at import time, and THIS file imports nothing at all.

WHAT THIS PACKAGE IS
====================
The scientific instrument in this repository is 372 Python files, but the *production
import closure* traced by AST from `s9/final.py` (the final held-out evaluation) is 18
modules -- 28 once VQE/CVaR are included.  Everything else is experiment debris.  `core`
is that closure, consolidated, made parallel, resumable and measurable, WITHOUT changing a
single scientific decision.

    core.pipeline   the four-stage pipeline (retrieve -> filter -> synthesise -> relax)
    core.bench      the benchmark harness that measures it, per stage, against the
                    `s9/final.py` reference path

    core.geometry   backbone build / Kabsch / the manifold projection
    core.data       the peptide database, folds, identity clusters
    core.predict    the distogram predictor
    core.energy     the classical energy terms
    core.amber      the ff14SB/GBn2 restrained relaxation
    core.quantum    the VQE / CVaR optimiser

THE ONE-LINE SWITCH
===================
`core.pipeline` never imports `protein_geometry` (or any other root module) directly.  It
asks `core.backend("geometry")`, and this module hands back `core.geometry` if it exists
AND exposes the whole contract, otherwise the root module it replaces.  So a consolidated
module going live is *literally* the act of committing the file -- no edit here, no edit
in `pipeline.py`.

The contract check is deliberate.  A partially-written `core/geometry.py` appearing on
disk mid-sprint must not silently break the pipeline or, far worse, silently change a
number; `backend()` requires every symbol in `CONTRACT[name]` to be present before it will
prefer the consolidated module, and records what it chose in `ACTIVE` so `core.bench` can
write the live backend set into every results file AND into every cache key.  A run whose
backends differ is a run with a different cache key: an optimised result can never be
served out of a baseline cache.

Set `CORE_BACKENDS=legacy` in the environment to force every fallback -- this is how the
`--baseline` arm of the harness guarantees it is measuring the reference path -- or
`CORE_BACKENDS=geometry,amber` to force a specific subset consolidated.
"""
from __future__ import annotations

import importlib
import os

__all__ = ["backend", "backend_name", "ACTIVE", "CONTRACT", "REPLACES", "backend_report"]

#: name -> (consolidated module, the root module it replaces)
REPLACES = {
    "geometry": ("core.geometry", "protein_geometry"),
    "data": ("core.data", "peptide_db"),
    "predict": ("core.predict", "distogram"),
    "energy": ("core.energy", "energy_terms"),
    #: the Legacy field's 11 decomposed terms and their fitted weights.  `core.energy`
    #: consolidates `energy_terms` AND `legacy_field`, so this entry's legacy fallback is
    #: the second of those, not the first.
    "legacy": ("core.energy", "legacy_field"),
    "amber": ("core.amber", "amber_refine"),
    "quantum": ("core.quantum", "foldvqe"),
    #: the numerics `s7/audit.py` owns (BLOSUM encode, batched Kabsch, pair distances).
    #: These are geometry, so a consolidated `core.geometry` may claim them; until it does
    #: they come from `s7.audit`, which is inside the production closure.
    "numerics": ("core.geometry", "s7.audit"),
    #: STAGE 3b, the manifold projection, consolidated in its own module: it carries the
    #: fold-disciplined Ramachandran prior as well as the optimiser, and neither belongs in
    #: `core.geometry`.  `s8.project` remains the reference the equivalence is measured
    #: against (`verify/project_equiv.json`).
    "project": ("core.project", "s8.project"),
}

#: What a consolidated module must expose before it is preferred over the module it
#: replaces.  Symbol-for-symbol identical call signatures are the switch's whole premise.
CONTRACT = {
    "geometry": ("build_backbone_batch", "extract_torsions", "kabsch_superpose", "rmsd"),
    "data": ("load", "folds", "clusters", "benchmark", "dev_set", "identity",
             "IDENTITY_THRESHOLD"),
    "predict": ("_fold_fragments", "train_fold", "Distogram", "_model_path"),
    "energy": (),
    "legacy": ("BatchLegacy", "TERMS", "FITTED_WEIGHTS", "SUBSETS"),
    "amber": ("refine_coords", "clear_cache"),
    "quantum": ("run_cvar_vqe", "cvar_exact", "free_energy"),
    "numerics": ("B62", "encode", "kabsch_rmsd_batch", "pair_index", "pair_dists"),
    "project": ("make_penalty", "lam_path"),
}

#: name -> the module string actually in use, filled in as backends are resolved.
ACTIVE: dict = {}

_CACHE: dict = {}


def _forced() -> set:
    v = os.environ.get("CORE_BACKENDS", "").strip().lower()
    if not v:
        return set(REPLACES)          # try everything
    if v in ("legacy", "none", "off", "baseline"):
        return set()
    return {x.strip() for x in v.split(",") if x.strip()}


def backend(name: str):
    """The live implementation of one backend, consolidated if it is ready, else legacy."""
    if name in _CACHE:
        return _CACHE[name]
    if name not in REPLACES:
        raise KeyError(f"unknown backend {name!r}; known: {sorted(REPLACES)}")
    opt, legacy = REPLACES[name]
    mod = None
    if name in _forced():
        try:
            cand = importlib.import_module(opt)
            missing = [s for s in CONTRACT[name] if not hasattr(cand, s)]
            if missing:
                if os.environ.get("CORE_VERBOSE"):
                    print(f"[core] {opt} present but incomplete "
                          f"(missing {missing}); using {legacy}", flush=True)
            else:
                mod = cand
                ACTIVE[name] = opt
        except ImportError:
            mod = None
    if mod is None:
        mod = importlib.import_module(legacy)
        ACTIVE[name] = legacy
    _CACHE[name] = mod
    return mod


def backend_name(name: str) -> str:
    backend(name)
    return ACTIVE[name]


def backend_report() -> dict:
    """`{backend: module actually used}` for every backend, for the results file."""
    return {n: backend_name(n) for n in sorted(REPLACES)}
