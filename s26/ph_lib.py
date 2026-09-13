"""s26/ph_lib.py -- shared helpers for the S26 Physics lane (PH).

Owned by lane PH.  Everything the three experiment scripts share lives here so that the
phase gate, the native-free loaders and the provenance stamp are one implementation.

THE PHASE GATE.  `s26/LANE_CONTRACT.md` section 4: no endpoint experiment (anything that
reads an RMSD to a native) may run until the coordinator posts `PHASE 0 SIGNED OFF` in
`s26/LEDGER.md`.  `require_gate()` is called by every function that reads a native
coordinate, a native torsion or a stored RMSD, and it raises if the line is absent.  A
census that reads native OMEGA angles only (ledger L5) goes through `native_backbone`, which
is labelled ORACLE DIAGNOSTIC and is the one native reader the coordinator allowed before
the gate.

NATIVE-FREE LOADERS.  `s8/generate_univ/<pdb>.npz` carries `rr` (ORACLE per-window RMSD)
and `nat_ca` (ORACLE native trace) beside the window bank.  `univ_nativefree` opens the
npz and reads ONLY the keys asked for, never those two.  `prod_record_nativefree` strips
every stored RMSD from the production record.  The oracle twins call `require_gate` first.

BASIS DISCIPLINE.  Point cloud (`rmsd_avg`, 3.0483), built chain (`rmsd_arm`, 3.2148),
relaxed chain (`rmsd_full`, 3.2355).  Every function that returns an RMSD names its basis
in its docstring; every artefact row carries a `basis` key.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys
import time
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)
LEDGER = os.path.join(HERE, "LEDGER.md")
GATE_LINE = "PHASE 0 SIGNED OFF"
UNIV = os.path.join(ROOT, "s8", "generate_univ")
CACHE_AMB = os.path.join(ROOT, "s24", "cache_amber")
PROD_KEY = "1fc9f2dcf489e2fb"
PROD_DIR = os.path.join(ROOT, "bench_results", "cache", PROD_KEY)
PDB_DIRS = (os.path.join(ROOT, "pdbs_ext"), os.path.join(ROOT, "pdbs"))
SALT = "s26ph"

#: keys of the production record that are RMSDs to the native or derived from one.
ORACLE_KEYS = ("shipped", "pool_best", "pool_mean", "top_m_best", "top_m_mean",
               "rmsd_avg", "rmsd_fit", "rmsd_arm", "rmsd_full", "d_rmsd",
               "rmsd_fit_single")
#: keys of the window universe that read the native.
ORACLE_UNIV_KEYS = ("rr", "nat_ca")

K = 500
M = 75
IDEAL_CA_CA = 3.803954938363982      # the built chain's virtual bond (e_reproduce.json)
CIS_OMEGA_DEG = 30.0                 # |omega| below this is cis (brief section 3.2)
CIS_CA_CA = 3.3                      # consecutive CA-CA below this is cis (trans 3.80, cis 2.9)


# ------------------------------------------------------------------ the phase gate
def gate_open(ledger: str = LEDGER) -> bool:
    """True iff the ledger carries the sign-off as a HEADING or as a standalone line.

    A substring test is not enough: ledger L5 quotes the phrase inside a sentence
    ('... before "PHASE 0 SIGNED OFF". The cis-peptide census ...'), and that must not open
    the gate.  Found by `grep -c` returning 1 before any sign-off existed.
    """
    try:
        with open(ledger, encoding="utf-8") as fh:
            for line in fh:
                s = line.strip()
                if s.startswith("#") and GATE_LINE in s:
                    return True
                if s.startswith(GATE_LINE):
                    return True
    except OSError:
        return False
    return False


def require_gate(what: str, ledger: str = LEDGER) -> None:
    """Raise unless the coordinator has posted the sign-off line."""
    if not gate_open(ledger):
        raise RuntimeError(
            f"PHASE GATE CLOSED: {what!r} reads a native RMSD and `{GATE_LINE}` is not in "
            f"{os.path.relpath(ledger, ROOT)}. Run the native-free mode instead.")


# ------------------------------------------------------------------ targets and folds
def targets() -> List[Dict]:
    """The 126 dev targets in pinned pdb order (pdb, n, fold, seq). Reads no native."""
    from s12 import instrument as I
    return I.targets()


def pdb_list() -> List[str]:
    return [t["pdb"] for t in targets()]


def folds_of(pdbs: Sequence[str]) -> np.ndarray:
    from s24 import stats_lib as ST
    return ST.pinned_folds(list(pdbs))


# ------------------------------------------------------------------ native-free loaders
def univ_nativefree(pdb: str, keys: Iterable[str] = ("W", "PHI", "PSI", "order", "sim",
                                                     "org", "S")) -> Dict:
    """Window universe with ONLY the requested keys read; `rr` and `nat_ca` are refused."""
    bad = [k for k in keys if k in ORACLE_UNIV_KEYS]
    if bad:
        raise ValueError(f"native-free loader asked for oracle keys {bad}")
    z = np.load(os.path.join(UNIV, f"{pdb}.npz"), allow_pickle=True)
    out = {"pdb": str(z["pdb"]), "seq": str(z["seq"]), "n": int(z["n"]),
           "fold": int(z["fold"])}
    for k in keys:
        v = z[k]
        if k in ("W", "PHI", "PSI"):
            v = np.asarray(v, np.float64)
        out[k] = v
    return out


def univ_oracle(pdb: str) -> Dict:
    """`s12.instrument.load_univ`, gated. ORACLE: carries `rr` and `nat_ca`."""
    require_gate(f"univ_oracle({pdb})")
    from s12 import instrument as I
    return I.load_univ(pdb)


def pool_idx_from_order(u: Dict, k: int = K) -> np.ndarray:
    """Indices into the universe of the shipped K=500 BLOSUM pool (`I.pool_idx`, no native)."""
    return np.asarray(u["order"][:k], int)


def prod_record_nativefree(pdb: str) -> Dict:
    """The production record with every stored RMSD removed."""
    with open(os.path.join(PROD_DIR, f"{pdb}.json"), encoding="utf-8") as fh:
        rec = json.load(fh)
    for k in ORACLE_KEYS:
        rec.pop(k, None)
    return rec


def prod_record_oracle(pdb: str) -> Dict:
    require_gate(f"prod_record_oracle({pdb})")
    with open(os.path.join(PROD_DIR, f"{pdb}.json"), encoding="utf-8") as fh:
        return json.load(fh)


def cache_amber(pdb: str) -> Dict:
    """`s24/cache_amber/<pdb>.npz`: 500 genuine ff14SB/GBn2 single points, native-free."""
    z = np.load(os.path.join(CACHE_AMB, f"{pdb}.npz"), allow_pickle=True)
    return {"pdb": str(z["pdb"]), "n": int(z["n"]), "fold": int(z["fold"]), "k": int(z["k"]),
            "e_amber": np.asarray(z["e_amber"], float),
            "score_dist": np.asarray(z["score_dist"], float),
            "universe_idx": np.asarray(z["universe_idx"], int),
            "amber_verify_max_rel": float(z["amber_verify_max_rel"])}


# ------------------------------------------------------------------ the native (ORACLE)
def native_pdb_path(pdb: str) -> str:
    """The deposited file of ONE dev target. Never enumerates a directory."""
    for d in PDB_DIRS:
        for name in (pdb + ".pdb", pdb.lower() + ".pdb", pdb.upper() + ".pdb"):
            p = os.path.join(d, name)
            if os.path.exists(p):
                return p
    raise FileNotFoundError(f"no deposited file for {pdb} under {PDB_DIRS}")


def native_backbone(pdb: str, model_index: int = 0) -> Tuple[str, np.ndarray, np.ndarray,
                                                            np.ndarray]:
    """ORACLE DIAGNOSTIC. `(seq, N, CA, C)` of one deposited model (default model 1, the
    model the instrument scores against). Allowed before the gate ONLY for the omega
    census (ledger L5); nothing here computes an RMSD."""
    from core import geometry as geo
    seq, N, CA, C = geo.parse_pdb(native_pdb_path(pdb), model_index=model_index)
    return str(seq), np.asarray(N, float), np.asarray(CA, float), np.asarray(C, float)


def native_ensemble(pdb: str) -> List[Tuple[str, np.ndarray, np.ndarray, np.ndarray]]:
    """ORACLE DIAGNOSTIC. Every deposited model (NMR ensembles)."""
    from core import geometry as geo
    return [(str(s), np.asarray(N, float), np.asarray(CA, float), np.asarray(C, float))
            for s, N, CA, C in geo.parse_pdb_ensemble(native_pdb_path(pdb))]


# ------------------------------------------------------------------ geometry helpers
def omega_deg(CA: np.ndarray, C: np.ndarray, N: np.ndarray) -> np.ndarray:
    """omega_i = dihedral(CA_i, C_i, N_{i+1}, CA_{i+1}) in degrees, i = 0..n-2."""
    from core import geometry as geo
    CA = np.asarray(CA, float); C = np.asarray(C, float); N = np.asarray(N, float)
    if len(CA) < 2:
        return np.zeros(0)
    return np.degrees(geo.dihedral_batch(CA[:-1], C[:-1], N[1:], CA[1:]))


def consecutive_ca(CA: np.ndarray) -> np.ndarray:
    CA = np.asarray(CA, float)
    return np.linalg.norm(CA[..., 1:, :] - CA[..., :-1, :], axis=-1)


def rg_of(X: np.ndarray) -> np.ndarray:
    """Radius of gyration of every (n,3) in a (b,n,3) stack (or one (n,3))."""
    X = np.asarray(X, float)
    single = X.ndim == 2
    if single:
        X = X[None]
    c = X - X.mean(1, keepdims=True)
    r = np.sqrt((c ** 2).sum(-1).mean(1))
    return r[0] if single else r


def min_sep3(X: np.ndarray) -> np.ndarray:
    """Minimum CA-CA distance over pairs with |i-j| >= 3, per structure."""
    X = np.asarray(X, float)
    single = X.ndim == 2
    if single:
        X = X[None]
    n = X.shape[1]
    i, j = np.triu_indices(n, k=3)
    if len(i) == 0:
        out = np.full(len(X), np.inf)
    else:
        d = np.linalg.norm(X[:, i, :] - X[:, j, :], axis=-1)
        out = d.min(1)
    return out[0] if single else out


def ca_contacts(X: np.ndarray, cutoff: float, min_sep: int = 3) -> np.ndarray:
    """Number of CA pairs closer than `cutoff` at |i-j| >= min_sep, per structure."""
    X = np.asarray(X, float)
    single = X.ndim == 2
    if single:
        X = X[None]
    n = X.shape[1]
    i, j = np.triu_indices(n, k=min_sep)
    if len(i) == 0:
        out = np.zeros(len(X), int)
    else:
        d = np.linalg.norm(X[:, i, :] - X[:, j, :], axis=-1)
        out = (d < cutoff).sum(1)
    return out[0] if single else out


def superpose_onto(X: np.ndarray, T: np.ndarray) -> np.ndarray:
    """X (n,3) superposed onto T (n,3), proper rotations only (frozen instrument)."""
    from s12 import instrument as I
    return I.superpose_batch(np.asarray(X, float)[None], np.asarray(T, float))[0]


def ca_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    from s12 import instrument as I
    return I.ca_rmsd(np.asarray(a, float), np.asarray(b, float))


def rigid_basis(CA: np.ndarray) -> np.ndarray:
    """Orthonormal (3n, 6) basis of rigid-body directions at CA (`s15.align_lib.rigid_basis`)."""
    from s15 import align_lib as AL
    return AL.rigid_basis(np.asarray(CA, float))


def random_displacement(CA: np.ndarray, magnitude: float, rng) -> np.ndarray:
    """S16's matched-magnitude control, reproduced: an isotropic Gaussian direction on the
    CA trace with the six rigid-body components removed, scaled so that
    ||g|| / sqrt(n) == magnitude (the per-atom RMS displacement). Native-free."""
    A = np.asarray(CA, float)
    n = len(A)
    g = rng.normal(size=(n, 3))
    Q = rigid_basis(A)
    g = (g.ravel() - Q @ (Q.T @ g.ravel())).reshape(n, 3)
    g *= float(magnitude) * math.sqrt(n) / max(np.linalg.norm(g), 1e-12)
    return g


def stable_rng(*parts):
    from s15 import seed as SD
    return SD.stable_rng(*parts, salt=SALT)


def argsort_stable(x) -> np.ndarray:
    return np.argsort(np.asarray(x, float), kind="stable")


# ------------------------------------------------------------------ provenance / io
def save(path: str, obj: Dict, rows: Optional[List[Dict]] = None,
         complete_keys: Optional[Sequence[str]] = None, n_expected: Optional[int] = None,
         module_file: Optional[str] = None) -> str:
    from s24 import stats_lib as ST
    return ST.save_atomic(path, obj, complete_keys=complete_keys, rows=rows,
                          n_expected=n_expected, module_file=module_file)


def cell_dir(name: str) -> str:
    d = os.path.join(RESULTS, f"{name}_cells")
    os.makedirs(d, exist_ok=True)
    return d


def write_cell(name: str, pdb: str, row: Dict) -> str:
    """Atomic per-target checkpoint so a governor kill loses nothing."""
    p = os.path.join(cell_dir(name), f"{pdb}.json")
    t = p + f".tmp{os.getpid()}"
    with open(t, "w", encoding="utf-8") as fh:
        json.dump(row, fh, default=_jsonable)
    os.replace(t, p)
    return p


def read_cells(name: str) -> Dict[str, Dict]:
    out = {}
    for p in sorted(glob.glob(os.path.join(cell_dir(name), "*.json"))):
        with open(p, encoding="utf-8") as fh:
            r = json.load(fh)
        out[r["pdb"]] = r
    return out


def _jsonable(o):
    if hasattr(o, "tolist"):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return str(o)


def mean_se(x) -> Dict[str, float]:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {"mean": float("nan"), "se": float("nan"), "n": 0}
    return {"mean": float(x.mean()), "se": float(x.std(ddof=1) / math.sqrt(len(x))) if len(x) > 1
            else float("nan"), "median": float(np.median(x)), "n": int(len(x))}


class Clock:
    def __init__(self):
        self.t0 = time.time()

    def __call__(self) -> float:
        return time.time() - self.t0
