"""One deterministic persistent cache for the data, geometry and prediction layers.

The rule this file exists to enforce: a cache key contains **every parameter that can
change the result**, including the version of the code that produced it. A scientifically
invalid collision -- two different computations reading each other's answer -- is far worse
than a slow run, and it is silent. So:

* The key is a BLAKE2b digest over a canonical encoding of ``(namespace, version, params)``.
  Canonical means sorted keys, explicit type tags, and numpy arrays hashed by dtype, shape
  and raw bytes -- ``1`` and ``1.0`` and ``"1"`` all key differently, and a float32 array
  never collides with the float64 array of the same values.
* Every entry is written with a **sidecar of the parameters that produced it**. On read the
  sidecar is compared against the requested parameters, and a mismatch raises
  `CacheCollision` instead of returning the wrong array. A digest collision is
  astronomically unlikely; a *bug in the key construction* -- a parameter someone forgot to
  pass -- is not, and that is what the sidecar actually catches.
* `version` is a per-namespace integer bumped by hand when the implementation changes in a
  way that changes results. Bumping it orphans the old entries rather than corrupting them.

THE 1.5 GB LANDMINE. ``esm_cache.npz`` is an object-array npz; ``np.load`` on it
materialises the whole 1.5 GB, and doing that inside a long-lived process has repeatedly
taken this box to 94-96% RAM. `extract_subset` runs the load in a **throwaway subprocess
that exits**, so the peak is that process's and it is returned to the OS immediately. Never
call `np.load` on a multi-gigabyte object npz from a process that has to keep running.

Storage layout::

    core_cache/<namespace>/<key[:2]>/<key>.npz     arrays
    core_cache/<namespace>/<key[:2]>/<key>.json    {"namespace","version","params","written"}

Writes are atomic (temp file + ``os.replace``), so an interrupted run leaves no half-file
for the next one to read.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.environ.get("CORE_CACHE", os.path.join(BASE, "core_cache"))

#: Set ``CORE_CACHE_OFF=1`` to bypass the disk layer entirely (tests that must not read a
#: stale entry, and any run whose whole point is to recompute).
DISABLED = os.environ.get("CORE_CACHE_OFF", "") == "1"

_MEM: Dict[str, Any] = {}
_MEM_ORDER: list = []
_MEM_MAX = int(os.environ.get("CORE_CACHE_MEM", "64"))

STATS: Dict[str, Dict[str, int]] = {}


class CacheCollision(RuntimeError):
    """A key resolved to an entry whose recorded parameters differ from the request."""


# ------------------------------------------------------------------ key construction
def _canon(obj: Any) -> Any:
    """Canonical, type-tagged, JSON-encodable form of a cache parameter.

    Type tags are what stop ``1``, ``1.0``, ``True`` and ``"1"`` from sharing a key: JSON
    alone renders the first three as ``1``/``1.0``/``true`` but a bare tuple and list
    identically, and ``np.int64(1)`` not at all.
    """
    if obj is None or isinstance(obj, bool):
        return ["b", obj]
    if isinstance(obj, (int, np.integer)):
        return ["i", int(obj)]
    if isinstance(obj, (float, np.floating)):
        # repr of a float is round-trip exact in Python 3, so this distinguishes values
        # that differ in the last bit -- which a tolerance-based comparison would not.
        return ["f", repr(float(obj))]
    if isinstance(obj, str):
        return ["s", obj]
    if isinstance(obj, bytes):
        return ["y", hashlib.blake2b(obj, digest_size=16).hexdigest()]
    if isinstance(obj, np.ndarray):
        a = np.ascontiguousarray(obj)
        return ["a", str(a.dtype), list(a.shape),
                hashlib.blake2b(a.tobytes(), digest_size=16).hexdigest()]
    if isinstance(obj, dict):
        return ["d", [[_canon(k), _canon(v)] for k, v in sorted(obj.items(),
                                                                key=lambda kv: str(kv[0]))]]
    if isinstance(obj, tuple):
        return ["t", [_canon(v) for v in obj]]
    if isinstance(obj, (list, set, frozenset)):
        if isinstance(obj, (set, frozenset)):
            return ["e", sorted(json.dumps(_canon(v), sort_keys=True) for v in obj)]
        return ["l", [_canon(v) for v in obj]]
    # Anything else is refused rather than str()-ed: str() of an object with a default
    # __repr__ embeds its memory address, which would make the key nondeterministic across
    # runs -- a cache that never hits is merely slow, but one that hits on an address reuse
    # is wrong.
    raise TypeError(f"cache parameter of type {type(obj).__name__} is not canonicalisable; "
                    f"pass a primitive, array, or container of them")


def key(namespace: str, version: int, **params: Any) -> str:
    """Deterministic hex key for ``(namespace, version, params)``."""
    blob = json.dumps(["k", namespace, int(version), _canon(params)],
                      sort_keys=True, separators=(",", ":"))
    return hashlib.blake2b(blob.encode("utf-8"), digest_size=20).hexdigest()


def _paths(namespace: str, k: str) -> Tuple[str, str]:
    d = os.path.join(ROOT, namespace, k[:2])
    return os.path.join(d, k + ".npz"), os.path.join(d, k + ".json")


def _bump(namespace: str, field: str) -> None:
    s = STATS.setdefault(namespace, {"hit": 0, "miss": 0, "mem_hit": 0, "write": 0})
    s[field] += 1


def stats() -> Dict[str, Dict[str, int]]:
    """Per-namespace ``{hit, miss, mem_hit, write}`` counters for this process."""
    return {k: dict(v) for k, v in STATS.items()}


def hit_rate() -> Dict[str, float]:
    out = {}
    for ns, s in STATS.items():
        tot = s["hit"] + s["mem_hit"] + s["miss"]
        out[ns] = (s["hit"] + s["mem_hit"]) / tot if tot else 0.0
    return out


def reset_stats() -> None:
    STATS.clear()


# ------------------------------------------------------------------ store
def _mem_put(k: str, value: Any) -> None:
    if k in _MEM:
        return
    _MEM[k] = value
    _MEM_ORDER.append(k)
    while len(_MEM_ORDER) > _MEM_MAX:
        _MEM.pop(_MEM_ORDER.pop(0), None)


def load(namespace: str, version: int, **params: Any) -> Optional[Dict[str, np.ndarray]]:
    """Return the cached arrays for these parameters, or None.

    Raises `CacheCollision` if an entry exists under this key whose recorded parameters
    differ from the ones asked for -- which means the key is missing a parameter.
    """
    k = key(namespace, version, **params)
    if k in _MEM:
        _bump(namespace, "mem_hit")
        return _MEM[k]
    if DISABLED:
        _bump(namespace, "miss")
        return None
    npz, side = _paths(namespace, k)
    if not (os.path.exists(npz) and os.path.exists(side)):
        _bump(namespace, "miss")
        return None
    try:
        with open(side) as fh:
            meta = json.load(fh)
    except Exception:                                          # pragma: no cover
        _bump(namespace, "miss")
        return None
    want = json.dumps(_canon(params), sort_keys=True, separators=(",", ":"))
    if meta.get("params") != want:
        raise CacheCollision(
            f"{namespace}/{k}: stored parameters differ from the request. The key is "
            f"missing a parameter that changes the result.\n  stored:    {meta.get('params')}"
            f"\n  requested: {want}")
    z = np.load(npz, allow_pickle=True)
    out = {kk: z[kk] for kk in z.files}
    z.close()
    _mem_put(k, out)
    _bump(namespace, "hit")
    return out


def store(namespace: str, version: int, arrays: Dict[str, np.ndarray],
          **params: Any) -> str:
    """Write ``arrays`` under these parameters. Atomic. Returns the key."""
    k = key(namespace, version, **params)
    _mem_put(k, {kk: np.asarray(v) for kk, v in arrays.items()})
    if DISABLED:
        return k
    npz, side = _paths(namespace, k)
    os.makedirs(os.path.dirname(npz), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(npz), suffix=".npz")
    os.close(fd)
    try:
        np.savez(tmp, **arrays)
        # np.savez appends .npz to a path that lacks it; mkstemp gave us one that has it.
        os.replace(tmp, npz)
    except Exception:                                          # pragma: no cover
        for p in (tmp, tmp + ".npz"):
            if os.path.exists(p):
                os.remove(p)
        raise
    meta = {"namespace": namespace, "version": int(version),
            "params": json.dumps(_canon(params), sort_keys=True, separators=(",", ":")),
            "written": time.time()}
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(side), suffix=".json")
    with os.fdopen(fd, "w") as fh:
        json.dump(meta, fh)
    os.replace(tmp, side)
    _bump(namespace, "write")
    return k


def cached(namespace: str, version: int, compute, **params: Any) -> Dict[str, np.ndarray]:
    """`load` if present, else ``compute()`` -> dict of arrays -> `store` -> return it."""
    got = load(namespace, version, **params)
    if got is not None:
        return got
    out = compute()
    if not isinstance(out, dict):
        raise TypeError("cached(compute=...) must return a dict of arrays")
    store(namespace, version, out, **params)
    return out


def clear(namespace: Optional[str] = None) -> int:
    """Delete a namespace (or the whole cache). Returns the number of entries removed."""
    _MEM.clear()
    _MEM_ORDER.clear()
    target = ROOT if namespace is None else os.path.join(ROOT, namespace)
    n = 0
    for root, _dirs, files in os.walk(target, topdown=False):
        for f in files:
            if f.endswith((".npz", ".json")) and f != ".gitignore":
                os.remove(os.path.join(root, f))
                n += f.endswith(".npz")
        if root != ROOT and not os.listdir(root):
            os.rmdir(root)
    return n


def size_on_disk() -> Dict[str, Tuple[int, int]]:
    """``{namespace: (entries, bytes)}``."""
    out: Dict[str, Tuple[int, int]] = {}
    if not os.path.isdir(ROOT):
        return out
    for ns in sorted(os.listdir(ROOT)):
        d = os.path.join(ROOT, ns)
        if not os.path.isdir(d):
            continue
        n = b = 0
        for root, _dirs, files in os.walk(d):
            for f in files:
                if f.endswith(".npz"):
                    n += 1
                    b += os.path.getsize(os.path.join(root, f))
        out[ns] = (n, b)
    return out

# ------------------------------------------------------------------ the 1.5 GB landmine
_EXTRACT = r"""
import json, sys
import numpy as np
src, dst, keyfile = sys.argv[1], sys.argv[2], sys.argv[3]
with open(keyfile) as fh:
    want = set(json.load(fh))
z = np.load(src, allow_pickle=True)
seqs = [str(s) for s in z["seqs"]]
vals = z["vals"]
idx = [i for i, s in enumerate(seqs) if s in want]
np.savez_compressed(dst, seqs=np.array([seqs[i] for i in idx], dtype=object),
                    vals=np.array([vals[i] for i in idx], dtype=object))
print(len(idx))
"""


def extract_subset(src_npz: str, sequences: Sequence[str], dst_npz: str,
                   timeout: int = 1800) -> int:
    """Copy just ``sequences`` out of a large object-array npz, in a throwaway subprocess.

    ``esm_cache.npz`` is ~1.5 GB and ``np.load(..., allow_pickle=True)`` materialises all
    of it. Doing that inside a process that then keeps running has repeatedly pushed this
    16 GB box to 94-96%. Here the load happens in a child that exits immediately
    afterwards, so the peak is transient and the parent never grows. Caching the resulting
    per-target table instead of the whole bank cut one downstream footprint from 2.0 GB to
    0.26 GB.

    Returns the number of sequences written.
    """
    if not os.path.exists(src_npz):
        raise FileNotFoundError(src_npz)
    fd, keyfile = tempfile.mkstemp(suffix=".json")
    with os.fdopen(fd, "w") as fh:
        json.dump(sorted(dict.fromkeys(sequences)), fh)
    fd, script = tempfile.mkstemp(suffix=".py")
    with os.fdopen(fd, "w") as fh:
        fh.write(_EXTRACT)
    try:
        os.makedirs(os.path.dirname(os.path.abspath(dst_npz)), exist_ok=True)
        r = subprocess.run([sys.executable, script, src_npz, dst_npz, keyfile],
                           capture_output=True, text=True, timeout=timeout)
        if r.returncode != 0:                                  # pragma: no cover
            raise RuntimeError(f"extract_subset failed: {r.stderr[-2000:]}")
        return int(r.stdout.strip().splitlines()[-1])
    finally:
        for p in (keyfile, script):
            if os.path.exists(p):
                os.remove(p)

# ------------------------------------------------------------------ box discipline
import ctypes                                                   # noqa: E402


class _MEMSTAT(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]


def mem_pct() -> Optional[int]:
    """Physical memory load via `GlobalMemoryStatusEx`, read in process.

    NOT ``Get-CimInstance Win32_OperatingSystem``, which costs minutes per call on a loaded
    box and turns the gate protecting the box into the slowest thing on it.
    """
    if os.name != "nt":                                        # pragma: no cover
        return None
    s = _MEMSTAT()
    s.dwLength = ctypes.sizeof(_MEMSTAT)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s)):   # pragma: no cover
        return None
    return int(s.dwMemoryLoad)


def mem_gate(tag: str = "", start: int = 86, stop: int = 92, wait: int = 25,
             tries: int = 12) -> Optional[int]:
    """Report above `start`, block above `stop`. Three sibling agents share this box."""
    for k in range(tries):
        m = mem_pct()
        if m is None or m < stop:
            if m is not None and m >= start:
                print(f"  [mem {tag}] {m}%", flush=True)
            return m
        print(f"  [mem gate {tag}] {m}% >= {stop}%, waiting {wait}s ({k + 1}/{tries})",
              flush=True)
        time.sleep(wait)
    return mem_pct()                                           # pragma: no cover
