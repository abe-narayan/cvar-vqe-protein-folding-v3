"""S25 / PHYSICS LANE -- shared machinery for the 7-configuration comparison suite.

Pre-registered in `s25/PREREG_PHYS.md`. Read that file before reading this one; the six
Rule-0 forks and the four falsifiers live there, not in the code.

WHAT THIS MODULE OWNS
=====================
* `zrank` and `zmoment` -- the two normalisations, and `combine()`, which is the ONLY place a
  configuration's energy is built. Every configuration goes through the same three lines, so
  "same normalisation across configurations" is a property of the code, not of the caller's
  discipline.
* `CONFIGS` -- the seven channel subsets, in the brief's order, frozen.
* `channels(pdb)` -- the three raw channels for one target, from ONE pool object, with the
  AMBER cache's pool identity and distogram score asserted bit-exact before they are used.
* `AmberLock` on `s25/results/LOCK_AMBER`, announced on open and release.

WHY THE ENERGY IS BUILT HERE AND NOWHERE ELSE
=============================================
`control-must-match-the-operators-space` is the most repeated error in this project's memory
(three instances in two sprints). The 7-configuration suite's whole value is that the seven
arms differ in EXACTLY one thing. So the pool, the readout, the selector settings and the
normalisation are constructed once, in this module, and the configuration is a set of channel
NAMES -- there is no API path by which one configuration can receive a different pool, a
different normalisation or a different readout from another.

NORMALISATION, RESTATED FROM THE PRE-REGISTRATION
=================================================
    zrank(x) = (rankdata(x) - mean) / sd              strictly monotone; mean 0, sd 1 exactly
    E_S      = zrank( sum_{c in S} zrank(x_c) )       equal unit weight on matched marginals

The OUTER zrank is not decoration. Without it a |S|-channel sum has sd
`sqrt(|S| + |S|(|S|-1) rho_bar)` -- 1.00 / 1.41 / 1.73 at rho_bar ~ 0 for |S| = 1/2/3 -- and
the CVaR objective trades energy against `T*H`, so the temperature would silently differ
between configurations. `zrank` is idempotent, so single-channel configurations are unchanged.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE_LEG = os.path.join(HERE, "cache_legacy")
CACHE_AMB = os.path.join(ROOT, "s24", "cache_amber")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE_LEG, exist_ok=True)

SALT = "s25phys"
K = 500
M_PROD = 75

#: pinned by the native-free tail-size probe in PREREG_PHYS 1.3 -- median realised |tail| = 75
#: at alpha = 0.18, which is the production rung. NEVER chosen on RMSD.
ALPHA = 0.18
TEMP = 0.5
LAYERS = 3
ITERS = 80
SEED = 0

CHANNELS = ("LEG", "AMB", "DIS")
#: the brief's order, frozen. `name` is what appears in every table.
CONFIGS: Tuple[Tuple[int, str, Tuple[str, ...]], ...] = (
    (1, "Legacy",                    ("LEG",)),
    (2, "AMBER",                     ("AMB",)),
    (3, "Distogram",                 ("DIS",)),
    (4, "Legacy+Distogram",          ("LEG", "DIS")),
    (5, "AMBER+Distogram",           ("AMB", "DIS")),
    (6, "Legacy+AMBER",              ("LEG", "AMB")),
    (7, "Legacy+AMBER+Distogram",    ("LEG", "AMB", "DIS")),
)


# ======================================================================= normalisation
def zrank(x) -> np.ndarray:
    """Standardised rank. Strictly monotone -- no ordering, argmin or level set changes."""
    from scipy.stats import rankdata
    r = rankdata(np.asarray(x, float))
    return (r - r.mean()) / max(float(r.std()), 1e-12)


def zmoment(x) -> np.ndarray:
    """Raw moment z. THE DECLARED SECONDARY. Clash-dominated on AMBER by construction."""
    v = np.asarray(x, float)
    return (v - v.mean()) / max(float(v.std()), 1e-12)


def combine(ch: Dict[str, np.ndarray], subset, norm: str = "rank") -> np.ndarray:
    """THE ONLY place a configuration's energy is built. Equal unit weight on matched marginals.

    `norm="rank"`  -> E_S = zrank( sum zrank(x_c) )      the PRIMARY
    `norm="moment"`-> E_S = sum zmoment(x_c)             the DECLARED SECONDARY

    The secondary deliberately does NOT re-standardise: its whole purpose is to show what the
    unconditioned scale does, and hiding that behind an outer standardisation would defeat it.
    """
    if norm == "rank":
        s = np.zeros(len(ch[subset[0]]), float)
        for c in subset:
            s = s + zrank(ch[c])
        return zrank(s)
    if norm == "moment":
        s = np.zeros(len(ch[subset[0]]), float)
        for c in subset:
            s = s + zmoment(ch[c])
        return s
    raise ValueError(f"unknown normalisation {norm!r}")


# =========================================================================== the channels
def legacy_cached(cand) -> np.ndarray:
    """Genuine 11-term Legacy at DEFAULT_WEIGHTS, cached per target. Never fitted."""
    f = os.path.join(CACHE_LEG, f"{cand.pdb}.npz")
    if os.path.exists(f):
        z = np.load(f)
        if int(z["k"]) == cand.k and np.array_equal(np.asarray(z["universe_idx"], int),
                                                    np.asarray(cand.meta["universe_idx"], int)):
            return np.asarray(z["e_leg"], float)
    e = H.score_legacy(cand)
    tmp = f + f".tmp{os.getpid()}"
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, pdb=cand.pdb, k=cand.k, e_leg=np.asarray(e, float),
                            universe_idx=np.asarray(cand.meta["universe_idx"], int))
    os.replace(tmp, f)
    return np.asarray(e, float)


def channels(pdb: str) -> Tuple[object, Dict[str, np.ndarray]]:
    """The one pool object and its three raw channels, with the cache identity ASSERTED.

    Two assertions, both cheap and both catching a class of error this project has already
    paid for: the AMBER cache must have been computed on THIS pool (`universe_idx` bit-equal)
    and the cached distogram score must reproduce bit-for-bit under the current code. A stale
    cache that silently disagrees with the live pool would corrupt every AMBER configuration.
    """
    z = np.load(os.path.join(CACHE_AMB, f"{pdb}.npz"), allow_pickle=True)
    cand = H.Candidates.from_universe(pdb, k=int(z["k"]))
    ui = np.asarray(z["universe_idx"], int)
    assert np.array_equal(np.asarray(cand.meta["universe_idx"], int), ui), \
        f"{pdb}: AMBER cache was computed on a DIFFERENT pool -- refusing to use it"
    sc = H.score_shipped(cand)
    assert np.array_equal(sc, np.asarray(z["score_dist"], float)), \
        f"{pdb}: cached distogram score does not reproduce bit-for-bit -- cache is stale"
    assert float(z["amber_verify_max_rel"]) == 0.0, \
        f"{pdb}: cached amber_verify_max_rel != 0"
    return cand, {"LEG": legacy_cached(cand), "AMB": np.asarray(z["e_amber"], float),
                  "DIS": sc}


# ============================================================================== readout
def rmsd_of_set(cand, idx) -> float:
    """THE readout, used by every arm and every configuration. POINT-CLOUD basis.

    ORACLE evaluation (`nat_ca`) of an ACHIEVABLE, native-free selection.
    """
    C, _ = I.coordinate_average(cand.W[np.asarray(idx, int)])
    return float(I.ca_rmsd(C, cand.nat_ca))


# ============================================================================== locking
class AmberLock:
    """Exclusive `s25/results/LOCK_AMBER`. Announced on open and on release, bounded holds."""

    def __init__(self, tag: str = "", wait: float = 10.0, tries: int = 180):
        self.path = os.path.join(RESULTS, "LOCK_AMBER")
        self.tag, self.wait, self.tries, self.held = tag, float(wait), int(tries), False

    def __enter__(self):
        for _ in range(self.tries):
            try:
                fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, json.dumps(dict(pid=os.getpid(), tag=self.tag, lane="PHYS",
                                             opened=time.strftime("%Y-%m-%d %H:%M:%S"))).encode())
                os.close(fd)
                self.held = True
                print(f"[LOCK_AMBER] OPENED by lane PHYS pid={os.getpid()} tag={self.tag}",
                      flush=True)
                return self
            except FileExistsError:
                time.sleep(self.wait)
        raise TimeoutError("LOCK_AMBER held by another process")

    def __exit__(self, *exc):
        if self.held:
            try:
                os.remove(self.path)
            except OSError:
                pass
            print(f"[LOCK_AMBER] RELEASED by lane PHYS pid={os.getpid()} tag={self.tag}",
                  flush=True)
        return False


def targets() -> List[str]:
    return sorted(t["pdb"] for t in I.targets())
