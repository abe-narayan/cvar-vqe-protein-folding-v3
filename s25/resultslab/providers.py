"""s25/resultslab/providers.py -- WHERE PREDICTIONS COME FROM.  Architecture-agnostic.

A provider is a callable `f(pdb) -> dict | None`.

    None                       the configuration has NO result for that target.  ABSENT.
                               It is written as absent.  It is never filled with an estimate.
    {"ca": (n,3) array, ...}   a prediction, plus whatever provenance the producer knows.

Recognised provenance keys (all optional; anything unknown is carried through verbatim into
the record's `extra` block):

    selector          e.g. "cvar_vqe", "classical_argmin"     -- WHO chose
    hamiltonians      list, e.g. ["legacy", "amber"]          -- WHAT scored
    distogram_used    bool
    cvar_alpha        float
    qubits            int
    ansatz_layers     int
    candidate_count   int      candidates the selector saw
    retained_count    int      candidates that survived into the emitted average
    seed              int

WHY A REGISTRY AND NOT A SWITCH.  The Sprint-25 architecture is not frozen.  Nothing in this
package may know which score, which selector or which Hamiltonian combination wins.  A
research lane freezes, drops a file, registers a name, and the whole lab -- structures,
schema, leaderboard, renderer -- follows without an edit here.  If you find yourself adding
an `if configuration == ...` to this module, you have re-introduced the coupling.

THE TWO BUILT-IN PROVIDERS ARE BOTH CONTROLS, NOT RESULTS.
  * `incumbent()` replays the SHIPPED point cloud (`avg_ca`, mean 3.0483 A over the 126 dev
    targets).  It is the baseline every configuration is paired against.  Its selector is
    the shipped distogram Bayes-risk argmin -- a CLASSICAL ranker, so it is labelled as one.
  * `synthetic()` is native-plus-seeded-noise.  It exists ONLY to prove the pipeline end to
    end before real results exist.  It is marked `is_synthetic_provider`, and `register`
    REFUSES to bind it to any real configuration name -- it can only be attached to a
    `_`-prefixed test fixture, whose structures are quarantined under
    `results/structures/_fixtures/`.

PROVENANCE DEFAULTS TO UNKNOWN, NOT TO GENUINE.  `from_json`/`from_npz` used to stamp a
loaded structure `is_synthetic: False` when the file said nothing at all.  A guard whose
default is "this is real" is not a guard.  Absent provenance is now `UNKNOWN` and
`schema.build(require_provenance=True)` -- which the frozen build uses -- refuses it.

THE GATE THAT DOES NOT DEPEND ON ANY OF THIS.  `exportlib.pool_oracle_gate` is the release
condition: no configuration may score below the ORACLE best member of its own K=500
candidate pool (mean 1.7108 A).  It is native-free, machine-checkable and needs nobody to
remember a flag, and it catches every failure the two rules above catch plus the ones
nobody thought of.
"""
from __future__ import annotations

import json
import os
import sys
import zlib
from typing import Callable, Dict, List, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                     # noqa: E402
from s25.resultslab import exportlib as EX          # noqa: E402

#: Sentinel for "this producer did not say".  Distinct from False (declared genuine) and
#: True (declared synthetic).  `schema.build(require_provenance=True)`, which the frozen
#: build uses, refuses to write a record whose provenance is UNKNOWN.
UNKNOWN = "unknown"

Provider = Callable[[str], Optional[Dict[str, object]]]

_REGISTRY: Dict[str, Provider] = {}
_META: Dict[str, Dict[str, object]] = {}


def register(configuration: str, provider: Provider, **static_meta) -> None:
    """Bind a configuration LABEL to a provider.  `static_meta` is merged into every record.

    TWO STRUCTURAL REFUSALS, both of them lessons rather than taste:

    * a label must be one of the pinned configuration names or a `_`-prefixed TEST FIXTURE;
    * a SYNTHETIC provider may never be bound to a real configuration name.  The proof build
      used to register native-plus-noise under six of the seven real names, which meant the
      only thing standing between a synthetic structure and a page reading "legacy: 2.10 A"
      was a status flag someone had to remember to look at.  It is now impossible to express.
    """
    if configuration not in EX.CONFIGURATIONS and configuration not in EX.FIXTURES:
        raise KeyError("configuration %r is not one of the pinned labels %s (fixtures: %s)"
                       % (configuration, EX.CONFIGURATIONS, EX.FIXTURES))
    if getattr(provider, "is_synthetic_provider", False) and configuration in EX.CONFIGURATIONS:
        raise ValueError(
            "REFUSED: a synthetic provider cannot be registered under the real configuration "
            "name %r. Synthetic arms belong under a test-fixture label (%s), which is "
            "quarantined in results/structures/_fixtures/ and can never be mistaken for a "
            "result." % (configuration, ", ".join(EX.FIXTURES)))
    _REGISTRY[configuration] = provider
    _META[configuration] = dict(static_meta)


def registered() -> List[str]:
    """Registered labels: pinned configurations in display order, then any test fixtures."""
    return ([c for c in EX.CONFIGURATIONS if c in _REGISTRY]
            + [c for c in EX.FIXTURES if c in _REGISTRY])


def get(configuration: str) -> Provider:
    return _REGISTRY[configuration]


def static_meta(configuration: str) -> Dict[str, object]:
    return dict(_META.get(configuration, {}))


def clear() -> None:
    _REGISTRY.clear()
    _META.clear()


# ------------------------------------------------------------------ built-in: incumbent
def incumbent(key: str = I.PROD_KEY, field: str = "avg_ca") -> Provider:
    """Replay a cached production run.  `field` selects the BASIS-CARRYING array.

    `avg_ca` is the point cloud (3.0483 A).  `fit_ca` is the lam=0 built chain (3.2041 A)
    and `ca` the lam=0.3 ramah arm (3.2148 A) -- but those are CHAIN coordinates and must be
    exported on the built-chain basis, not relabelled as point clouds.
    """
    def f(pdb: str):
        p = os.path.join(ROOT, "bench_results", "cache", key, "%s.json" % str(pdb).upper())
        if not os.path.exists(p):
            return None
        with open(p) as fh:
            r = json.load(fh)
        if r.get(field) is None:
            return None
        return {"ca": np.asarray(r[field], float),
                "selector": "distogram_bayes_risk_argmin (classical control)",
                "hamiltonians": [],
                "distogram_used": True,
                "cvar_alpha": None,
                "qubits": None,
                "ansatz_layers": None,
                "candidate_count": int(r.get("n_candidates") or 0) or None,
                "retained_count": int(r.get("n_top") or 0) or None,
                "seed": None,
                "source": "bench_results/cache/%s/%s.json:%s" % (key, str(pdb).upper(), field),
                "is_synthetic": False}
    return f


# ------------------------------------------------------------------ built-in: synthetic
def synthetic(sigma: float = 1.6, seed: int = 0, tag: str = "") -> Provider:
    """Native + seeded isotropic displacement.  A PIPELINE PROOF, never a result.

    The displacement is drawn from a per-target stable RNG so the same call reproduces the
    same structures on any machine.
    """
    def f(pdb: str):
        nat = EX.native_ca(pdb)
        # NOT `hash()`: Python salts string hashing per process (PYTHONHASHSEED), so a
        # hash-seeded RNG makes a "deterministic" build produce different structures on
        # every run.  Caught by `test_synthetic_provider_is_deterministic`.
        key = "%s|%s|%d" % (str(pdb).upper(), tag, int(seed))
        rng = np.random.default_rng(zlib.crc32(key.encode()) & 0xFFFFFFFF)
        # Deliberately emits ONLY what it knows.  Everything else is left to the static
        # metadata the caller registered, so a synthetic arm still displays the labels of
        # the configuration it is standing in for instead of blanking them.
        return {"ca": nat + rng.normal(0.0, float(sigma), nat.shape),
                "selector": "SYNTHETIC (no selector)",
                "seed": int(seed),
                "source": "providers.synthetic(sigma=%.3f, seed=%d, tag=%r)" % (sigma, seed, tag),
                "is_synthetic": True}
    f.is_synthetic_provider = True     # `register` refuses to bind this to a real name
    return f


# ------------------------------------------------------------------ built-in: from a file
def from_json(path: str, **static) -> Provider:
    """The hand-off a frozen research lane uses.  No code change here, just a file.

    Accepted shapes, both keyed by PDB id:
        {"1A13": [[x,y,z], ...], ...}
        {"1A13": {"ca": [[x,y,z], ...], "cvar_alpha": 0.15, ...}, ...}
    A target absent from the file is ABSENT, not zero and not interpolated.
    """
    with open(path) as fh:
        blob = json.load(fh)
    blob = blob.get("targets", blob) if isinstance(blob, dict) else blob

    def f(pdb: str):
        got = blob.get(str(pdb).upper())
        if got is None:
            return None
        if isinstance(got, dict):
            if got.get("ca") is None:
                return None
            out = {k: v for k, v in got.items() if k != "ca"}
            out["ca"] = np.asarray(got["ca"], float)
        else:
            out = {"ca": np.asarray(got, float)}
        out.setdefault("source", os.path.relpath(path, ROOT).replace("\\", "/"))
        for k, v in static.items():
            out.setdefault(k, v)
        # PROVENANCE IS NOT ASSUMED GENUINE.  This used to default `is_synthetic` to False,
        # so a loaded structure carrying no provenance at all was silently stamped REAL.
        # For a leakage guard the safe default is the opposite: absent provenance is
        # UNKNOWN, and `schema.build(require_provenance=True)` refuses to ship it.
        out.setdefault("is_synthetic", UNKNOWN)
        return out
    return f


def from_npz(path: str, key_ca: str = "ca", key_pdb: str = "pdb", **static) -> Provider:
    """As `from_json`, for an npz holding parallel `pdb` (str) and `ca` ((T, n, 3) object) arrays."""
    z = np.load(path, allow_pickle=True)
    pdbs = [str(p).upper() for p in z[key_pdb]]
    cas = list(z[key_ca])
    blob = {p: np.asarray(c, float) for p, c in zip(pdbs, cas)}

    def f(pdb: str):
        got = blob.get(str(pdb).upper())
        if got is None:
            return None
        out = {"ca": got, "source": os.path.relpath(path, ROOT).replace("\\", "/")}
        out.update(static)
        out.setdefault("is_synthetic", UNKNOWN)   # absent provenance is UNKNOWN, not genuine
        return out
    return f
