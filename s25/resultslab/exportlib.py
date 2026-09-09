"""s25/resultslab/exportlib.py -- THE STRUCTURE EXPORT LAYER.  Deterministic, reusable.

WHAT THIS IS.  One place that turns `(target, configuration, CA trace)` into a valid,
self-describing PDB file, and one place that reads it back.  Everything downstream of the
results lab -- the schema, the leaderboard, the renderer -- reads files written here, so
the single property that makes the whole lab trustworthy is asserted here and nowhere else:

    THE RMSD RECOMPUTED FROM AN EXPORTED FILE REPRODUCES THE INSTRUMENT'S OWN
    `s12.instrument.ca_rmsd` TO WITHIN PDB'S %8.3f QUANTISATION.

`verify_roundtrip()` is that assertion.  If it fails, nothing downstream means anything.

-------------------------------------------------------------------------------------
THE BASIS.  RULED, s25 LEDGER L4/L8.  READ THIS BEFORE QUOTING ANY NUMBER FROM HERE.
-------------------------------------------------------------------------------------
Measured over the 126 pinned dev targets from `bench_results/cache/1fc9f2dcf489e2fb`:

    POINT CLOUD   `avg_ca`   the coordinate average, no ideal-geometry chain   3.0483 A
    BUILT CHAIN   `fit_ca`   that average projected to the nearest ideal chain 3.2041 A
    BUILT CHAIN   `ca`       ... with the ramah prior at lam=0.3  <- WE SHIP THIS  3.2148 A

**THE BUILT CHAIN IS THE PRODUCTION RESULT, at lam=0.3 -- SYNTHESIS stage 3, the arm this
project ships, incumbent 3.2148 A.  The point cloud is retained everywhere as an explicitly
labelled NON-PHYSICAL INTERMEDIATE, and the two are always reported PAIRED.**

The projection gap for the SYNTHESIS arm is +0.1664 A (median +0.0977, IQR [+0.005, +0.338]),
and the built chain is actually BETTER on 16 of 126 targets.  The often-quoted +0.156 belongs
to the lam=0 fit arm, which is a DIFFERENT arm and is not what we ship.

The point cloud is not a molecule.  Its mean virtual Ca-Ca bond is 2.9614 A against a
native 3.8122 A -- 22.3% contracted, worst single bond 0.649 A, shorter than a covalent
C-C bond, with 80 of 126 targets under 3.4 A.  `ca_bond_stats` carries that measurement on
every exported file and every record, so the label is self-evidencing.

This module's own round-trip gate is what made the ruling implementable rather than merely
arguable: **the point-cloud arm cannot pass the gate, because a 0.649 A "bond" is not a
structure that can be written as a valid PDB and round-tripped as one.**  There is no file
to verify.  `export_pair` therefore emits both bases from one call, so there is no way to
obtain one number here without also obtaining the other.

The brief's Phase-II note quotes a built-chain figure of 3.2126 A.  **That number does not
reproduce anywhere in this repository**; the reproducible figures are 3.2041 (lam=0) and
3.2148 (lam=0.3).  Nothing here quotes 3.2126.

Structural separation, so a basis mix-up is visible rather than silent:
    results/structures/chain/T001__legacy.pdb   built chain  (PRODUCTION)
    results/structures/T001__legacy.pdb         point cloud  (intermediate)
    results/structures/T001__native.pdb         the shared Ca referent
    <system temp>/s25_resultslab_fixtures/...   TEST FIXTURES ONLY, never under results/

-------------------------------------------------------------------------------------
RULE 0 -- THE SIX OPERATOR FORKS, EACH NAMING THE ALTERNATIVE NOT TAKEN
-------------------------------------------------------------------------------------
1. FUNCTIONAL.  Exported geometry is the CA trace the configuration emitted, verbatim.
   NOT TAKEN: any post-hoc superposition-onto-native baked into the file.  The file holds
   the prediction in its own frame; superposition happens in the viewer, at draw time, so
   the stored coordinates can never be silently native-informed.
2. BASIS.  Primary export is the BUILT CHAIN; the point cloud is exported beside it and
   labelled a non-physical intermediate (s25 L4/L8).  NOT TAKEN: keeping the point cloud as
   the headline, which reads 0.166 A better and cannot be written as a valid molecule.  The
   direction of the correction is the thing that makes it credible.
3. READOUT.  RMSD stored in a file is recomputed FROM THE QUANTISED COORDINATES that file
   contains.  NOT TAKEN: writing the in-memory float64 RMSD into the header, which would
   have made the header disagree with the file by up to ~1e-3 A and made round-trip
   verification vacuous.
4. NORMALISATION.  T-numbers are assigned from `s12.instrument.targets()` order, which is
   the pinned pdb-sorted glob of `s8/generate_univ/*.npz`.  NOT TAKEN: sorting by length,
   by fold, or by RMSD -- any of which would renumber targets whenever a result changes.
   The map is persisted with a sha256 over the ordered id list and refuses to move.
5. NULL.  The pipeline is proven end-to-end against a SYNTHETIC provider (native plus a
   seeded displacement) in a SANDBOX, under `_`-prefixed fixture labels, writing to the
   system temp directory.  NOT TAKEN: the earlier proof build, which registered synthetic
   arms under six of the seven REAL configuration names and wrote 1,134 PDBs into
   `results/structures/` whose headers said nothing about it.  The 1,134 files were removed in
   the 2026-09-09 clean; the post-mortem that makes them citable is
   `docs/quarantine-synthetic-proofbuild.md`.
6. THE LABEL.  A configuration's name in this package is exactly the string the caller
   registered; nothing here decides what "production" means.  NOT TAKEN: hard-coding the
   current best pipeline as `production`, which is the mistake the sequencing rule for
   this sprint explicitly forbids.

-------------------------------------------------------------------------------------
REUSE
-------------------------------------------------------------------------------------
* `core.geometry.write_pdb` already emits correctly-columned ATOM records.  This module
  does NOT copy it: `_atom_lines` is asserted BYTE-IDENTICAL to it by
  `test_export.test_atom_lines_match_core_writer`.  A local writer exists only because
  `core.geometry.write_pdb` takes a single `remark` string and the lab needs a parseable
  REMARK 999 metadata block.
* `core.geometry`'s reader is NOT reusable for CA-only files -- `_read_models` drops any
  residue missing N, CA or C by design -- so `read_pdb` here is a minimal column reader.
  It is asserted against `core.geometry.parse_pdb` on backbone files.
* `s12.instrument.ca_rmsd` is the ONLY RMSD in this package.  Nothing is re-derived.
* `s24.stats_lib.save_atomic` stamps every JSON artefact.

-------------------------------------------------------------------------------------
THE FOUR RELEASE GATES.  ALL FOUR MUST PASS BEFORE ANY REAL EXPORT.
-------------------------------------------------------------------------------------
(a) `pool_oracle_gate`     mean RMSD >= the ORACLE best of its own K=500 pool (1.7108 A).
                           Native-free and a hard physical impossibility -- but the WEAKEST
                           of the four: it fires on little short of an exact copy, and every
                           synthetic arm in the quarantined proof build cleared it.
(b) `difficulty_gate`      corr(per-target RMSD, that target's pool_best) >= 0.30 AND
                           sd(per-target RMSD) >= 0.70.  THE ONE THAT DISCRIMINATES.
                           Measured here: incumbent +0.6988 / 1.6466 on the point cloud and
                           +0.7157 / 1.7361 on the production chain; synthetic arms
                           -0.12 .. +0.08 / 0.20 .. 0.47.  A real method is harder on hard
                           targets; isotropic noise about the native is not.
(c) `verify_header`        every PDB carries PROVENANCE / MODULE_SHA / GIT_COMMIT / BASIS in
                           its own REMARK 999 block, and the frozen build accepts only
                           `genuine`.  THE PRIMARY GATE: a label that does not travel with
                           the file is not a label.
(d) fixtures in temp       `FIXTURE_DIR` is the system temp directory and fixture labels are
                           `_`-prefixed; `providers.register` refuses to bind a synthetic
                           provider to a real configuration name at all.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import sys
import tempfile
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                     # noqa: E402
from s24 import stats_lib as ST                     # noqa: E402

# --------------------------------------------------------------------------- layout
RESULTS_DIR = os.path.join(ROOT, "results")
STRUCT_DIR = os.path.join(RESULTS_DIR, "structures")
CHAIN_DIR = os.path.join(STRUCT_DIR, "chain")
SUMMARY_DIR = os.path.join(RESULTS_DIR, "summary")
SITE_DIR = os.path.join(RESULTS_DIR, "site")
TARGET_MAP = os.path.join(SUMMARY_DIR, "target_map.json")

for _d in (STRUCT_DIR, CHAIN_DIR, SUMMARY_DIR, SITE_DIR):
    os.makedirs(_d, exist_ok=True)

# --------------------------------------------------------------------------- constants
#: THE PRODUCTION BASIS (s25 L4/L8 ruling).  The emitted point cloud is projected to the
#: nearest ideal-geometry backbone and THAT is the reported result.  It is a real chain: it
#: can be written as a valid PDB and it round-trips.
BASIS_BUILT_CHAIN = "built_chain_bb"

#: The NON-PHYSICAL INTERMEDIATE, retained and reported alongside, always labelled.  Its
#: virtual Ca-Ca bonds are 22.3% contracted (2.9614 A against a native 3.8122 A, worst single
#: bond 0.649 A), so it is a point cloud and not a molecule.
BASIS_POINT_CLOUD = "point_cloud_ca"

#: The referent both bases are scored against: `nat_ca` from the pinned universe.
BASIS_NATIVE = "native_ca"

PRIMARY_BASIS = BASIS_BUILT_CHAIN
SECONDARY_BASIS = BASIS_POINT_CLOUD

#: The seven pre-registered comparison configurations plus the production arm.  These are
#: LABELS ONLY -- this module attaches no meaning to them.  A provider must be registered
#: for a name before anything can be exported under it.
CONFIGURATIONS = (
    "legacy",
    "amber",
    "distogram",
    "legacy_distogram",
    "amber_distogram",
    "legacy_amber",
    "legacy_amber_distogram",
    "production",
)
NATIVE = "native"

#: TEST-ONLY labels.  Every one begins with "_", none collides with a real configuration
#: name, and `struct_path` routes all of them into `results/structures/_fixtures/`.
#:
#: WHY THIS EXISTS.  The round-trip fixture used to name its controls "distogram" (the
#: native itself, RMSD exactly 0.000000) and "legacy" (native + 1.5 A noise), so the test's
#: own printed table contained the line `T001  distogram  14  0.000000  0.000000`.  The test
#: was correct; the LABEL was a loaded gun.  One careless paste of that line is the most
#: damaging sentence this project could emit.  Controls now carry names that cannot be
#: mistaken for a result, and they cannot be written into the real structure directory.
FIXTURES = ("_zero", "_noise", "_incumbent")

#: Fixture structures are written to the SYSTEM TEMP DIRECTORY, never under `results/`.
#: `results/structures/_fixtures/` was not far enough: a directory that lives beside the
#: real exports is one `cp -r` from being shipped, and the quarantined proof build proved
#: that a file separated from its JSON reads as a genuine result.
FIXTURE_DIR = os.path.join(tempfile.gettempdir(), "s25_resultslab_fixtures")

#: The production chain arm.  `lam=0.3` with the `ramah` torsion prior is SYNTHESIS stage 3,
#: the arm this project ships (incumbent 3.2148 A); it reproduces the cached `ca` field
#: bit-exactly.  `lam=0.0` is the unconstrained fit arm (3.2041 A) and is NOT what we ship.
PRODUCTION_LAM = 0.3

#: `s12.instrument.selfcheck` pins this: the mean, over the 126 dev targets, of the ORACLE
#: best Ca-RMSD available anywhere in each target's own K=500 candidate pool.
POOL_BEST_MEAN = 1.7108

#: PDB stores coordinates as %8.3f.  Worst-case per-axis error is 5e-4 A, so the worst-case
#: per-atom displacement is sqrt(3)*5e-4 = 8.66e-4 A and RMSD -- an RMS over atoms of a
#: quantity bounded by that -- cannot exceed it either.  Kabsch is 1-Lipschitz in the point
#: set, so the superposed RMSD is bounded the same way.  Budget 2e-3 A, ~2.3x the bound.
RMSD_TOL = 2e-3

#: Provenance vocabulary carried in EVERY structure header.  UNKNOWN is the default for a
#: producer that did not declare, and it FAILS the frozen build.
PROV_GENUINE, PROV_SYNTHETIC, PROV_UNKNOWN = "genuine", "synthetic", "unknown"

#: The REMARK 999 keys every exported structure must carry.  `verify_header` checks them and
#: `export_prediction` cannot omit them.
#:
#: WHY THIS IS THE PRIMARY GATE.  The quarantined proof build wrote 1,134 PDB files, six of
#: seven configurations native-plus-noise, and NOT ONE carried a synthetic marker in its
#: header -- `is_synthetic` reached the JSON and never the artefact.  Separated from
#: results.json, `T001__distogram.pdb` reads as a genuine result at 0.000000 A.  A label that
#: does not travel with the file is not a label.
REQUIRED_HEADER_KEYS = ("provenance", "module_sha", "git_commit", "basis",
                        "target_id", "pdb_id", "configuration", "rmsd_ca", "referent")

_ONE_TO_THREE = None
_PROV_CACHE = None


def _prov() -> Dict[str, str]:
    """`module_sha` + `git_commit` for the header.  Computed once per process."""
    global _PROV_CACHE
    if _PROV_CACHE is None:
        got = ST.provenance(__file__)
        _PROV_CACHE = {"module_sha": str(got.get("source_sha256")),
                       "git_commit": str(got.get("git_commit"))}
    return _PROV_CACHE


def _one_to_three():
    global _ONE_TO_THREE
    if _ONE_TO_THREE is None:
        from core import geometry as geo
        _ONE_TO_THREE = dict(geo.ONE_TO_THREE)
    return _ONE_TO_THREE


# =========================================================================== T-numbers
def _ordered_pdbs() -> List[str]:
    """The pinned dev-target order.  `benchmark-and-folds-must-be-pinned`: read, never recompute."""
    return [t["pdb"] for t in I.targets()]


def _map_digest(pdbs: Sequence[str]) -> str:
    return hashlib.sha256("\n".join(pdbs).encode()).hexdigest()[:16]


def build_target_map(force: bool = False) -> Dict[str, object]:
    """Create or load the persisted PDB-id <-> T-number mapping.

    T-numbers are assigned 1-based over `s12.instrument.targets()` order, which is the glob
    of `s8/generate_univ/*.npz` sorted by filename, i.e. pdb-id-sorted.  That order is a
    property of the pinned instrument, not of any result, so it does not move when a
    configuration changes.

    The persisted file carries a sha256 of the ordered id list.  If the instrument's order
    ever changes, this function RAISES rather than silently renumbering -- every exported
    filename, every schema record and every rendered page keys off these numbers.
    """
    pdbs = _ordered_pdbs()
    digest = _map_digest(pdbs)
    if os.path.exists(TARGET_MAP) and not force:
        with open(TARGET_MAP) as fh:
            got = json.load(fh)
        if got.get("digest") != digest:
            raise RuntimeError(
                "TARGET MAP DRIFT: %s was built over a different target order "
                "(stored digest %s, instrument now %s).  Refusing to renumber: every "
                "exported artefact keys off these T-numbers.  Investigate the instrument "
                "before passing force=True." % (TARGET_MAP, got.get("digest"), digest))
        return got

    t_of_pdb, pdb_of_t, rows = {}, {}, []
    for k, t in enumerate(I.targets()):
        tid = "T%03d" % (k + 1)
        t_of_pdb[t["pdb"]] = tid
        pdb_of_t[tid] = t["pdb"]
        rows.append({"target_id": tid, "pdb_id": t["pdb"], "sequence": t["seq"],
                     "length": int(t["n"]), "fold": int(t["fold"])})
    obj = {"digest": digest, "n": len(rows), "order_source":
           "s12.instrument.targets() -- sorted glob of s8/generate_univ/*.npz",
           "t_of_pdb": t_of_pdb, "pdb_of_t": pdb_of_t, "targets": rows}
    ST.save_atomic(TARGET_MAP, obj, complete_keys=("target_id", "pdb_id", "sequence",
                                                   "length", "fold"),
                   rows=rows, n_expected=len(rows), module_file=__file__)
    with open(TARGET_MAP) as fh:
        return json.load(fh)


_MAP_CACHE: Optional[Dict[str, object]] = None


def target_map() -> Dict[str, object]:
    global _MAP_CACHE
    if _MAP_CACHE is None:
        _MAP_CACHE = build_target_map()
    return _MAP_CACHE


def tid_of(pdb: str) -> str:
    return target_map()["t_of_pdb"][str(pdb).upper()]


def pdb_of(tid: str) -> str:
    return target_map()["pdb_of_t"][str(tid).upper()]


def target_rows() -> List[Dict[str, object]]:
    return list(target_map()["targets"])


_ROW_CACHE: Optional[Dict[str, Dict[str, object]]] = None


def row_of(pdb: str) -> Dict[str, object]:
    """One target's pinned row, by PDB id.  Cached; the export path asks for it per structure."""
    global _ROW_CACHE
    if _ROW_CACHE is None:
        _ROW_CACHE = {r["pdb_id"]: r for r in target_rows()}
    return _ROW_CACHE[str(pdb).upper()]


def is_fixture(configuration: str) -> bool:
    return str(configuration).startswith("_")


def stem(pdb_or_tid: str, configuration: str) -> str:
    """`T001__legacy` -- the deterministic artefact stem.  Accepts either id form."""
    s = str(pdb_or_tid).upper()
    tid = s if re.fullmatch(r"T\d{3}", s) else tid_of(s)
    if configuration != NATIVE and configuration not in CONFIGURATIONS \
            and configuration not in FIXTURES:
        raise KeyError("unknown configuration %r; known: %s (test fixtures: %s)"
                       % (configuration, ", ".join(CONFIGURATIONS + (NATIVE,)),
                          ", ".join(FIXTURES)))
    return "%s__%s" % (tid, configuration)


def struct_path(pdb_or_tid: str, configuration: str,
                basis: str = BASIS_BUILT_CHAIN) -> str:
    """Where a structure lives.  Fixtures are quarantined; the two bases never share a dir."""
    name = stem(pdb_or_tid, configuration) + ".pdb"
    if is_fixture(configuration):
        return os.path.join(FIXTURE_DIR, basis, name)
    return os.path.join(CHAIN_DIR if basis == BASIS_BUILT_CHAIN else STRUCT_DIR, name)


# =========================================================================== PDB writing
def _atom_lines(sequence: str, coords: Dict[str, np.ndarray]) -> List[str]:
    """ATOM records, byte-identical to `core.geometry.write_pdb`'s body.

    Asserted, not assumed -- see `test_export.test_atom_lines_match_core_writer`.  Atom
    order per residue is N, CA, C, O; absent names are skipped, so passing only {"CA": ...}
    yields a CA-only trace.
    """
    o2t = _one_to_three()
    lines: List[str] = []
    serial = 1
    for i, aa in enumerate(sequence):
        res = o2t.get(aa, "GLY")
        for name, el in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
            if name not in coords:
                continue
            x, y, z = (float(v) for v in np.asarray(coords[name])[i])
            nm = (" " + name) if len(name) < 4 else name
            lines.append(f"ATOM  {serial:>5d} {nm:<4s} {res:>3s} A{i + 1:>4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
            serial += 1
    return lines


_META_KEY = re.compile(r"^REMARK 999 ([A-Z0-9_]+)\s(.*)$")


def _meta_lines(meta: Dict[str, object]) -> List[str]:
    out = []
    for k in sorted(meta):
        v = meta[k]
        if v is None:
            v = ""
        s = str(v).replace("\n", " ")
        out.append("REMARK 999 %s %s\n" % (str(k).upper(), s))
    return out


def quantise(sequence: str, coords: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """The coordinates a PDB file will actually contain, obtained by formatting and reparsing.

    NOT `np.round(x, 3)`: numpy rounds half-to-even on the binary value while `%8.3f` rounds
    the decimal representation, and the two disagree on ties.  Formatting and reading back is
    exact by construction, which is what lets `export_prediction` write the file ONCE.
    """
    lines = _atom_lines(sequence, coords)
    names = [k for k in ("N", "CA", "C", "O") if k in coords]
    out = {k: np.empty((len(sequence), 3)) for k in names}
    per = len(names)
    for i in range(len(sequence)):
        for j, k in enumerate(names):
            ln = lines[i * per + j]
            out[k][i] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return out


def write_structure(path: str, sequence: str, coords: Dict[str, np.ndarray],
                    meta: Dict[str, object]) -> str:
    """Write one structure file with a machine-readable `REMARK 999` metadata block.

    ONE write, tmp + os.replace.  `export_prediction` used to write twice -- once to obtain
    the quantised coordinates and once to stamp the RMSD computed from them -- so an
    interruption between the two left a file on disk whose header said `RMSD_CA PENDING`.
    One did (`T006__amber.pdb`).  `quantise` removes the need for the first write entirely.
    """
    n = len(sequence)
    for k, v in coords.items():
        a = np.asarray(v, float)
        if a.shape != (n, 3):
            raise ValueError("coords[%r] has shape %s, expected (%d, 3)" % (k, a.shape, n))
        if not np.isfinite(a).all():
            raise ValueError("coords[%r] is not finite" % k)
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    body = _atom_lines(sequence, coords)
    tmp = path + ".tmp.%d" % os.getpid()
    with open(tmp, "w", newline="\n") as fh:
        fh.writelines(_meta_lines(meta))
        fh.writelines(body)
        fh.write("TER\nEND\n")
    _replace_retry(tmp, path)
    return path


def verify_header(path: str, require: Sequence[str] = REQUIRED_HEADER_KEYS,
                  allow_provenance: Sequence[str] = (PROV_GENUINE,)) -> Dict[str, str]:
    """GATE (c).  Refuse a structure file that does not carry its own provenance.

    Reads the file back off disk -- not the dict that was meant to be written -- and checks
    that every required REMARK 999 key is present and that `PROVENANCE` is one the caller
    accepts.  The frozen build accepts only `genuine`; `unknown` and `synthetic` both fail.
    """
    got = read_pdb(path)["meta"]
    missing = [k for k in require if not str(got.get(k, "")).strip()]
    if missing:
        raise ProvenanceViolation(
            "PROVENANCE GATE FAILED: %s carries no %s in its REMARK 999 header. A label that "
            "does not travel with the file is not a label -- 1,134 quarantined proof-build "
            "PDBs each read as a genuine result for exactly this reason."
            % (os.path.basename(path), ", ".join(missing)))
    if got["provenance"] not in allow_provenance:
        raise ProvenanceViolation(
            "PROVENANCE GATE FAILED: %s is marked PROVENANCE=%s; this build accepts only %s."
            % (os.path.basename(path), got["provenance"], "/".join(allow_provenance)))
    return got


class ProvenanceViolation(AssertionError):
    """Raised when a structure file does not carry acceptable provenance in its own header."""


def _replace_retry(tmp: str, path: str, tries: int = 12, wait: float = 0.05) -> None:
    """`os.replace` with a bounded retry.  Windows only, and not paranoia.

    A full build writes ~2,000 small files in a couple of minutes; the on-access virus
    scanner and the search indexer both open a file the instant it is created, and while
    they hold it `os.replace` fails with WinError 5 (Access denied).  Observed here on
    T006__amber.pdb during a 1,008-structure run.  Retrying is correct because the write
    itself succeeded -- only the rename lost a race -- and a bounded loop still fails loudly
    if the cause is a genuine permission problem rather than a transient handle.
    """
    import time
    for k in range(tries):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            if k == tries - 1:
                raise
            time.sleep(wait * (k + 1))


# =========================================================================== PDB reading
def read_pdb(path: str) -> Dict[str, object]:
    """Minimal column reader: returns sequence, per-atom-name (n,3) arrays, and REMARK 999.

    `core.geometry`'s reader cannot be reused here: it drops any residue that is missing
    N, CA or C, which is every residue of a CA-only trace.  This reader is asserted against
    `core.geometry.parse_pdb` on backbone files by the test module.
    """
    from core import geometry as geo
    t2o = {v: k for k, v in geo.ONE_TO_THREE.items()}
    meta: Dict[str, str] = {}
    order: List[int] = []
    res: Dict[int, str] = {}
    atoms: Dict[str, Dict[int, Tuple[float, float, float]]] = {}
    with open(path) as fh:
        for raw in fh:
            m = _META_KEY.match(raw.rstrip("\n"))
            if m:
                meta[m.group(1).lower()] = m.group(2).strip()
                continue
            if raw[:6] != "ATOM  ":
                continue
            name = raw[12:16].strip()
            rn = raw[17:20].strip().upper()
            ri = int(raw[22:26])
            if ri not in res:
                res[ri] = rn
                order.append(ri)
            atoms.setdefault(name, {})[ri] = (float(raw[30:38]), float(raw[38:46]),
                                              float(raw[46:54]))
    seq = "".join(t2o.get(res[i], "X") for i in order)
    out: Dict[str, object] = {"sequence": seq, "meta": meta, "n": len(order),
                              "resseq": order, "path": os.path.abspath(path)}
    for name, d in atoms.items():
        out[name] = np.array([d[i] for i in order if i in d], float)
    return out


# =========================================================================== natives
_NATIVE_CACHE: Dict[str, np.ndarray] = {}


def native_ca(pdb: str) -> np.ndarray:
    """The evaluation referent: `nat_ca` from the pinned universe.  Nothing else is native here.

    This is deliberately NOT read from `pdbs/*.pdb`: the instrument scores against
    `u["nat_ca"]` (model 1 of the deposited ensemble, as cached at universe-build time), so
    the exported native must be that array and no other, or the displayed RMSD would not be
    the reported RMSD.

    Cached: `load_univ` reads a whole window universe (`W` is ~(10000, n, 3)) to get one
    (n, 3) array, and a full 8-configuration build asks for it ~3,000 times.  Returned
    read-only so a caller cannot mutate the referent for everyone else.
    """
    key = str(pdb).upper()
    got = _NATIVE_CACHE.get(key)
    if got is None:
        got = np.asarray(I.load_univ(key)["nat_ca"], float)
        got.setflags(write=False)
        _NATIVE_CACHE[key] = got
    return got


# =========================================================================== the pool oracle
POOL_BEST_FILE = os.path.join(SUMMARY_DIR, "pool_best.json")
_POOL_BEST: Optional[Dict[str, float]] = None


def pool_best(force: bool = False) -> Dict[str, float]:
    """ORACLE: per target, the best Ca-RMSD available anywhere in its own K=500 pool.

    Read from `u["rr"][order[:500]].min()` -- the same array `s12.instrument.selfcheck`
    averages to 1.7108 A.  This is a LABEL and is used for exactly one purpose: the release
    gate below.  It never enters an inference-time decision.

    Cached to `results/summary/pool_best.json` because computing it costs one full universe
    load per target.
    """
    global _POOL_BEST
    if _POOL_BEST is not None and not force:
        return _POOL_BEST
    if os.path.exists(POOL_BEST_FILE) and not force:
        with open(POOL_BEST_FILE) as fh:
            got = json.load(fh)
        if got.get("digest") == target_map()["digest"]:
            _POOL_BEST = {k: float(v) for k, v in got["pool_best"].items()}
            return _POOL_BEST
    out = {}
    for t in target_rows():
        u = I.load_univ(t["pdb_id"])
        out[t["pdb_id"]] = float(np.asarray(u["rr"], float)[I.pool_idx(u)].min())
    rows = [{"pdb_id": k, "pool_best": v} for k, v in sorted(out.items())]
    ST.save_atomic(POOL_BEST_FILE,
                   {"digest": target_map()["digest"], "pool_best": out,
                    "mean": float(np.mean(list(out.values()))),
                    "source": "s12.instrument: u['rr'][pool_idx(u)].min(), K=500",
                    "note": "ORACLE. Used only by the release gate; never at inference."},
                   complete_keys=("pdb_id", "pool_best"), rows=rows,
                   n_expected=len(target_rows()), module_file=__file__)
    _POOL_BEST = out
    return out


class PoolOracleViolation(AssertionError):
    """Raised when a configuration scores below the best structure in its own candidate pool."""


def pool_oracle_gate(per_target: Dict[str, float], configuration: str,
                     hard: bool = True) -> Dict[str, object]:
    """THE RELEASE GATE.  Native-free, machine-checkable, nobody has to remember a flag.

    No real method can beat the ORACLE best member of its own K=500 candidate pool, whose
    mean over the 126 dev targets is 1.7108 A.  A configuration that does has native
    information in it by some route -- a leaked label, a fixture wired to a real
    configuration name, a provider pointed at `nat_ca`.  The gate does not care which.

    AGGREGATE is the release condition and is HARD: mean(config) < 1.7108 fails the build.

    PER TARGET is reported but is only a WARNING, and the distinction is not squeamishness.
    The emitted structure is a coordinate AVERAGE over retained candidates and is therefore
    not itself a pool member, so it CAN legitimately land closer to the native than any
    single window on an individual target.  Averaging cannot do that systematically across
    126 targets, which is why the aggregate is the condition that binds.
    """
    pb = pool_best()
    names = [p for p in per_target if p in pb]
    if not names:
        return {"configuration": configuration, "checked": 0, "status": "not checked"}
    vals = np.array([per_target[p] for p in names], float)
    ref = np.array([pb[p] for p in names], float)
    viol = [(p, float(per_target[p]), float(pb[p])) for p in names if per_target[p] < pb[p]]
    out = {"configuration": configuration, "checked": int(len(names)),
           "mean": float(vals.mean()), "pool_best_mean": float(ref.mean()),
           "pool_best_mean_pinned": POOL_BEST_MEAN,
           "margin": float(vals.mean() - ref.mean()),
           "n_per_target_violations": len(viol),
           "per_target_violations": viol[:20],
           "aggregate_violation": bool(vals.mean() < ref.mean())}
    out["status"] = "FAIL" if out["aggregate_violation"] else (
        "WARN" if viol else "PASS")
    if hard and out["aggregate_violation"]:
        raise PoolOracleViolation(
            "POOL-ORACLE GATE FAILED for configuration %r: mean %.4f A is BELOW the mean "
            "ORACLE best of its own K=500 candidate pools (%.4f A over %d targets). No real "
            "method can do that. Something in this configuration has native information: a "
            "leaked label, a fixture registered under a real configuration name, or a "
            "provider reading nat_ca. %d individual targets also beat their own pool best."
            % (configuration, out["mean"], out["pool_best_mean"], out["checked"], len(viol)))
    return out


# =========================================================================== gate (b)
#: Floors for the DIFFICULTY-TRACKING gate.  Calibrated on the incumbent against every
#: synthetic arm the proof build produced; see `difficulty_gate`.
DIFFICULTY_CORR_FLOOR = 0.30
DIFFICULTY_SD_FLOOR = 0.70


class DifficultyGateViolation(AssertionError):
    """Raised when a configuration's errors do not track target difficulty."""


def difficulty_gate(per_target: Dict[str, float], configuration: str,
                    corr_floor: float = DIFFICULTY_CORR_FLOOR,
                    sd_floor: float = DIFFICULTY_SD_FLOOR,
                    hard: bool = True) -> Dict[str, object]:
    """GATE (b).  A real method must be HARDER ON HARD TARGETS.  Isotropic noise is not.

    Two statistics, both free from what the lab already computes:

        corr(per-target RMSD, that target's ORACLE pool_best)   >= corr_floor
        sd(per-target RMSD)                                     >= sd_floor

    THE CALIBRATION, measured on this instrument rather than assumed:

        incumbent (genuine)          corr +0.6988   sd 1.6400
        every synthetic arm          corr -0.08 .. +0.15   sd 0.26 .. 0.58

    The separation is not marginal, which is why fixed floors at 0.30 / 0.70 are defensible.
    The gate is nonetheless calibrated on ONE genuine arm, so state its limitation with it:
    a real method that were uniformly accurate across easy and hard targets would trip it.
    If a configuration fails here, that is a claim to investigate, not a verdict to hide --
    but it does not ship until someone has looked.

    This is the gate that actually discriminates.  The pool-oracle gate (a) is kept because
    it is a different kind of check -- a hard physical impossibility rather than a
    statistical signature -- but on its own it fires only on something close to an exact
    copy, and every synthetic arm in the quarantined proof build cleared it.
    """
    pb = pool_best()
    names = sorted(p for p in per_target if p in pb)
    if len(names) < 12:
        return {"configuration": configuration, "checked": len(names),
                "status": "not checked (fewer than 12 targets)"}
    v = np.array([per_target[p] for p in names], float)
    d = np.array([pb[p] for p in names], float)
    sd = float(v.std(ddof=1))
    corr = float(np.corrcoef(v, d)[0, 1]) if sd > 0 and d.std() > 0 else 0.0
    out = {"configuration": configuration, "checked": len(names),
           "corr_with_pool_best": corr, "sd": sd,
           "corr_floor": float(corr_floor), "sd_floor": float(sd_floor),
           "corr_ok": bool(corr >= corr_floor), "sd_ok": bool(sd >= sd_floor)}
    out["status"] = "PASS" if (out["corr_ok"] and out["sd_ok"]) else "FAIL"
    if hard and out["status"] == "FAIL":
        raise DifficultyGateViolation(
            "DIFFICULTY GATE FAILED for configuration %r over %d targets: "
            "corr(RMSD, pool_best) = %+.4f (floor %+.2f), sd = %.4f (floor %.2f). A genuine "
            "method is harder on hard targets and spreads accordingly; the incumbent sits at "
            "corr +0.6988 / sd 1.6400 while isotropic noise about the native sits at "
            "-0.08..+0.15 / 0.26..0.58. Either this configuration is not reading the target, "
            "or it is genuinely uniform across difficulty -- and the second needs a human."
            % (configuration, out["checked"], corr, corr_floor, sd, sd_floor))
    return out


# =========================================================================== geometry check
def ca_bond_stats(ca: np.ndarray) -> Dict[str, float]:
    """Virtual Ca-Ca bond lengths.  The evidence that a point cloud is not a chain.

    s25 L4/L8: the 3.0483 A point-cloud incumbent has a mean virtual Ca-Ca bond of 2.9614 A
    against a native 3.8122 A -- 22.3% contracted, worst single bond 0.649 A, shorter than a
    covalent C-C bond.  Carried on every record so the "non-physical intermediate" label is
    self-evidencing rather than an assertion the reader has to take on trust.
    """
    d = np.linalg.norm(np.diff(np.asarray(ca, float), axis=0), axis=1)
    if len(d) == 0:
        return {"ca_bond_mean": float("nan"), "ca_bond_min": float("nan"),
                "ca_bond_max": float("nan")}
    return {"ca_bond_mean": float(d.mean()), "ca_bond_min": float(d.min()),
            "ca_bond_max": float(d.max())}


# =========================================================================== chain build
def chain_from_ca(ca: np.ndarray, seq: str, fold: int, lam: float = PRODUCTION_LAM,
                  maxiter: int = 300, multi: bool = True) -> Dict[str, object]:
    """Point cloud -> nearest ideal-geometry backbone (production stage 3b).

    `lam=0.0` is the unconstrained fit arm -- the one whose incumbent mean is 3.2041 A.
    `lam=0.3` adds the `ramah` torsion prior and is the 3.2148 A arm.  Returns the CA trace,
    the torsions, and the full N/CA/C/O/CB backbone from `core.geometry.build_backbone`.
    """
    from core import project as pj
    from core import geometry as geo
    C = np.asarray(ca, float)
    pen = pj.make_penalty("ramah", seq, int(fold)) if lam > 0 else None
    lams = (0.0,) if lam == 0.0 else (0.0, float(lam))
    path = pj.lam_path(C, pen, lams, maxiter=maxiter, multi=multi, grad="exact")
    arm = path[float(lam)]
    phi = np.asarray(arm[1], float)
    psi = np.asarray(arm[2], float)
    bb = geo.build_backbone(phi, psi)
    return {"ca": np.asarray(arm[0], float), "phi": phi, "psi": psi,
            "backbone": {k: np.asarray(v, float) for k, v in bb.items()},
            "lam": float(lam)}


# =========================================================================== export API
def native_path(pdb: str, sandbox: bool = False) -> str:
    """Where the shared Ca referent lives.  `sandbox=True` keeps a selftest out of results/."""
    return os.path.join(FIXTURE_DIR if sandbox else STRUCT_DIR, stem(pdb, NATIVE) + ".pdb")


def export_native(pdb: str, sandbox: bool = False, **extra) -> Dict[str, object]:
    """Write `T###__native.pdb`: the Ca referent BOTH bases are scored against.

    One file, not one per basis.  The built chain and the point cloud are both compared to
    the same `nat_ca`, so a second copy could only ever be a way for them to disagree.
    """
    t = row_of(pdb)
    nat = native_ca(pdb)
    meta = {"target_id": t["target_id"], "pdb_id": t["pdb_id"], "configuration": NATIVE,
            "sequence": t["sequence"], "length": t["length"], "fold": t["fold"],
            "basis": BASIS_NATIVE, "basis_role": "REFERENT -- the deposited Ca trace",
            "rmsd_ca": "0.000000", "provenance": PROV_GENUINE,
            "referent": "s8/generate_univ nat_ca", "source": "s12.instrument.load_univ"}
    meta.update(_prov())
    meta.update(ca_bond_stats(nat))
    meta.update({k.lower(): v for k, v in extra.items()})
    path = native_path(pdb, sandbox=sandbox)
    write_structure(path, t["sequence"], {"CA": nat}, meta)
    return {"path": path, "meta": meta}


def export_prediction(pdb: str, configuration: str, ca: np.ndarray,
                      basis: str = PRIMARY_BASIS,
                      backbone: Optional[Dict[str, np.ndarray]] = None,
                      provenance: str = PROV_UNKNOWN,
                      **extra) -> Dict[str, object]:
    """Write one prediction file, ATOMICALLY, carrying its own provenance and its own RMSD.

    TWO PROPERTIES, both of them repairs of a real defect:

    ONE WRITE.  The header RMSD is computed from `quantise(...)` -- the coordinates the file
    will contain, obtained by formatting and reparsing -- so the file is written once,
    complete, via tmp + os.replace.  The previous version wrote twice, and an interruption
    between the two left `T006__amber.pdb` on disk with `RMSD_CA PENDING`.

    PROVENANCE TRAVELS WITH THE FILE.  `provenance` is a REQUIRED header field, defaulting to
    UNKNOWN.  The quarantined proof build wrote 1,134 PDBs of native-plus-noise and not one
    said so in its own header, because the flag lived only in the JSON.  Read on its own,
    `T001__distogram.pdb` announced a genuine 0.000000 A result.
    """
    t = row_of(pdb)
    ca = np.asarray(ca, float)
    if ca.shape != (t["length"], 3):
        raise ValueError("%s/%s: CA shape %s, expected (%d, 3)"
                         % (pdb, configuration, ca.shape, t["length"]))
    if not np.isfinite(ca).all():
        raise ValueError("%s/%s: CA trace is not finite -- refusing to export a structure "
                         "whose coordinates are not numbers" % (pdb, configuration))
    if provenance not in (PROV_GENUINE, PROV_SYNTHETIC, PROV_UNKNOWN):
        raise ValueError("provenance must be one of genuine/synthetic/unknown, got %r"
                         % provenance)
    coords = {"CA": ca}
    if basis == BASIS_BUILT_CHAIN:
        if backbone is None:
            raise ValueError("built-chain basis requires the full backbone")
        coords = {k: np.asarray(v, float) for k, v in backbone.items() if k in ("N", "CA", "C", "O")}
        if not np.allclose(coords["CA"], ca, atol=1e-9):
            raise ValueError("backbone CA disagrees with the CA trace being exported")

    # what the file WILL contain -- so the header can state its own RMSD before it is written
    qz = quantise(t["sequence"], coords)
    rmsd_file = I.ca_rmsd(qz["CA"], native_ca(pdb))

    path = struct_path(pdb, configuration, basis)
    meta = {"target_id": t["target_id"], "pdb_id": t["pdb_id"],
            "configuration": configuration, "sequence": t["sequence"],
            "length": t["length"], "fold": t["fold"], "basis": basis,
            "basis_role": ("PRODUCTION RESULT" if basis == PRIMARY_BASIS
                           else "NON-PHYSICAL INTERMEDIATE -- not a molecule"),
            "provenance": provenance,
            "referent": "s8/generate_univ nat_ca",
            "rmsd_ca": "%.6f" % rmsd_file,
            "rmsd_ca_memory": "%.6f" % I.ca_rmsd(ca, native_ca(pdb))}
    if is_fixture(configuration) or provenance == PROV_SYNTHETIC:
        tags = ([("TEST FIXTURE" if is_fixture(configuration) else None),
                 (provenance.upper() if provenance != PROV_GENUINE else None)])
        meta["warning"] = ("*** NOT A RESULT -- %s *** This file must never be read as a "
                           "prediction, quoted, or copied out of its directory."
                           % ", ".join(t for t in tags if t))
    meta.update(_prov())
    meta.update(ca_bond_stats(qz["CA"]))
    meta.update({k.lower(): v for k, v in extra.items()})
    write_structure(path, t["sequence"], coords, meta)
    return {"path": path, "meta": meta, "rmsd": rmsd_file,
            "rmsd_memory": float(meta["rmsd_ca_memory"]),
            "provenance": provenance,
            "basis": basis, "ca_bond": ca_bond_stats(qz["CA"])}


def export_pair(pdb: str, configuration: str, ca: np.ndarray, lam: float = PRODUCTION_LAM,
                provenance: str = PROV_UNKNOWN, **extra) -> Dict[str, object]:
    """Export BOTH bases from one emitted point cloud, and report them PAIRED.

    s25 L4/L8 ruling: the built chain is the production result; the point cloud is retained
    everywhere as an explicitly labelled non-physical intermediate.  Reporting them from a
    single call is what makes "always paired" a property of the code rather than a habit --
    there is no way to obtain one of these numbers here without also obtaining the other.
    """
    ca = np.asarray(ca, float)
    ch = chain_from_ca(ca, row_of(pdb)["sequence"], row_of(pdb)["fold"], lam=lam)
    prim = export_prediction(pdb, configuration, ch["ca"], basis=PRIMARY_BASIS,
                             backbone=ch["backbone"], provenance=provenance,
                             chain_lam="%.3f" % lam, **extra)
    sec = export_prediction(pdb, configuration, ca, basis=SECONDARY_BASIS,
                            provenance=provenance, **extra)
    return {"primary": prim, "secondary": sec, "chain": ch,
            "rmsd": prim["rmsd"], "basis": PRIMARY_BASIS,
            "rmsd_secondary": sec["rmsd"], "basis_secondary": SECONDARY_BASIS,
            "delta": float(prim["rmsd"] - sec["rmsd"])}


# =========================================================================== verification
def verify_roundtrip(pdb: str, configuration: str, ca: np.ndarray,
                     basis: str = PRIMARY_BASIS,
                     backbone: Optional[Dict[str, np.ndarray]] = None,
                     tol: float = RMSD_TOL) -> Dict[str, object]:
    """Export, re-read, recompute -- and compare to the instrument's own `ca_rmsd`.

    THE CHECK THAT MAKES THE LAB TRUSTWORTHY.  Returns every intermediate so a failure is
    diagnosable rather than merely loud.
    """
    ca = np.asarray(ca, float)
    nat = native_ca(pdb)
    truth = I.ca_rmsd(ca, nat)
    rec = export_prediction(pdb, configuration, ca, basis=basis, backbone=backbone)
    got = read_pdb(rec["path"])

    nat_rec = export_native(pdb, sandbox=is_fixture(configuration))
    nat_got = read_pdb(nat_rec["path"])
    from_files = I.ca_rmsd(got["CA"], nat_got["CA"])

    d_coord = float(np.abs(got["CA"] - ca).max())
    out = {"pdb": str(pdb).upper(), "target_id": tid_of(pdb),
           "configuration": configuration, "basis": basis,
           "n": int(len(ca)),
           "rmsd_instrument": float(truth),
           "rmsd_from_file_vs_memory_native": float(rec["rmsd"]),
           "rmsd_from_file_vs_file_native": float(from_files),
           "rmsd_in_header": float(got["meta"]["rmsd_ca"]),
           "abs_err_vs_instrument": abs(float(from_files) - float(truth)),
           "max_coord_abs_err": d_coord,
           "seq_roundtrip_ok": bool(got["sequence"] == got["meta"]["sequence"]),
           "header_matches_file": abs(float(got["meta"]["rmsd_ca"]) - float(rec["rmsd"])) < 1e-6,
           "tol": float(tol)}
    out.update(ca_bond_stats(got["CA"]))
    out["ok"] = bool(out["abs_err_vs_instrument"] < tol
                     and out["max_coord_abs_err"] <= 5.001e-4
                     and out["seq_roundtrip_ok"] and out["header_matches_file"])
    return out


def verify_all_built_chain(pdbs: Optional[Sequence[str]] = None,
                           lam: float = PRODUCTION_LAM,
                           configuration: str = "_incumbent",
                           verbose: bool = True) -> Dict[str, object]:
    """THE FULL-INSTRUMENT ROUND TRIP on the PRODUCTION basis.  All 126 targets, not a fixture.

    Runs the incumbent point cloud through the real stage-3b projection, exports the built
    chain, re-reads it, and asserts the RMSD off disk reproduces the instrument's.  Also
    records the virtual Ca-Ca bond statistics of both bases, which is the measurement that
    decided the basis ruling in the first place.

    Writes under the `_incumbent` FIXTURE label, so nothing here lands in the real structure
    directory or under a real configuration name.
    """
    rows = []
    tg = [t["pdb_id"] for t in target_rows()] if pdbs is None else [str(p).upper() for p in pdbs]
    for k, pdb in enumerate(tg):
        rec = I.shipped_record(pdb)
        pc = np.asarray(rec["avg_ca"], float)
        ch = chain_from_ca(pc, sequence_of(pdb), fold_of(pdb), lam=lam)
        r = verify_roundtrip(pdb, configuration, ch["ca"], basis=BASIS_BUILT_CHAIN,
                             backbone=ch["backbone"])
        r2 = verify_roundtrip(pdb, configuration, pc, basis=BASIS_POINT_CLOUD)
        r["point_cloud_rmsd"] = r2["rmsd_from_file_vs_file_native"]
        r["point_cloud_abs_err"] = r2["abs_err_vs_instrument"]
        r["point_cloud_ok"] = r2["ok"]
        r["point_cloud_bond_mean"] = r2["ca_bond_mean"]
        r["point_cloud_bond_min"] = r2["ca_bond_min"]
        r["native_bond_mean"] = ca_bond_stats(native_ca(pdb))["ca_bond_mean"]
        ref = rec["ca"] if abs(lam - PRODUCTION_LAM) < 1e-12 else rec["fit_ca"]
        r["vs_cached_arm"] = float(np.abs(ch["ca"] - np.asarray(ref, float)).max())
        rows.append(r)
        if verbose and (k + 1) % 25 == 0:
            print("    %d/%d" % (k + 1, len(tg)), flush=True)
    ok = all(r["ok"] and r["point_cloud_ok"] for r in rows)
    out = {
        "n": len(rows), "lam": float(lam), "all_ok": bool(ok),
        "worst_abs_err": max(r["abs_err_vs_instrument"] for r in rows),
        "worst_abs_err_point_cloud": max(r["point_cloud_abs_err"] for r in rows),
        "tolerance": RMSD_TOL,
        "mean_rmsd_built_chain": float(np.mean([r["rmsd_from_file_vs_file_native"] for r in rows])),
        "mean_rmsd_point_cloud": float(np.mean([r["point_cloud_rmsd"] for r in rows])),
        "worst_reprojection_vs_cached_arm": max(r["vs_cached_arm"] for r in rows),
        "mean_ca_bond_built_chain": float(np.mean([r["ca_bond_mean"] for r in rows])),
        "mean_ca_bond_point_cloud": float(np.mean([r["point_cloud_bond_mean"] for r in rows])),
        "mean_ca_bond_native": float(np.mean([r["native_bond_mean"] for r in rows])),
        "min_ca_bond_point_cloud": min(r["point_cloud_bond_min"] for r in rows),
        "min_ca_bond_built_chain": min(r["ca_bond_min"] for r in rows),
        "rows": rows,
    }
    out["basis_gap"] = out["mean_rmsd_built_chain"] - out["mean_rmsd_point_cloud"]
    return out


# =========================================================================== helpers
def sequence_of(pdb: str) -> str:
    return row_of(pdb)["sequence"]


def fold_of(pdb: str) -> int:
    return int(row_of(pdb)["fold"])


def length_of(pdb: str) -> int:
    return int(row_of(pdb)["length"])


def module_hash(*paths: str) -> str:
    """sha256 over the given source files, in the order given.  Goes into every record."""
    h = hashlib.sha256()
    for p in paths:
        with open(p, "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()[:16]


if __name__ == "__main__":
    m = build_target_map()
    print("target map: n=%d digest=%s -> %s" % (m["n"], m["digest"], TARGET_MAP))
    for r in m["targets"][:3] + m["targets"][-2:]:
        print("  %s  %s  n=%2d fold=%d  %s" % (r["target_id"], r["pdb_id"], r["length"],
                                               r["fold"], r["sequence"]))
