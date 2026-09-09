"""SPRINT 15 -- a STABLE seed, because `hash()` is not one.

The defect this fixes. Every multi-start module in this sprint seeded its RNG with

    rng = np.random.default_rng(hash(pdb) % (2 ** 32))

Python's string hash is **salted per process** unless `PYTHONHASHSEED` is set, and the sprint
environment lock records it as `null`. Two consecutive interpreters return
`hash('1A13') = 217586290314588545` and `2408026022170661001`. So every run drew a DIFFERENT set
of multi-start initialisations, and no multi-start number in this sprint was bit-reproducible on
re-running. Found by the RESTRAINT agent while checking that its ladder shared starts with the
coordinator's reference arm.

What it did and did not break, stated precisely because the distinction decides which results
survive:

  * **Within one process it is harmless.** Every arm in a single run shares one start set, so
    every paired comparison inside a single result file is start-matched and valid. All the
    sprint's within-file conclusions stand.
  * **Across processes it is not.** A figure from one module compared with a figure from another
    compares different start draws, and published constants are not reproducible on re-running.

`blake2b` is stable across processes, platforms and Python versions, which `hash()` is not.
"""
from __future__ import annotations

import hashlib


def stable_seed(*parts, salt: str = "s15") -> int:
    """A deterministic 32-bit seed from any hashable, printable parts.

    Order matters and is part of the key, so `stable_seed(pdb, "combined")` and
    `stable_seed("combined", pdb)` differ, exactly as the tuple hashes they replace did.
    """
    key = salt + "|" + "|".join(repr(p) for p in parts)
    return int.from_bytes(hashlib.blake2b(key.encode("utf-8"), digest_size=4).digest(), "big")


def stable_rng(*parts, salt: str = "s15"):
    import numpy as np
    return np.random.default_rng(stable_seed(*parts, salt=salt))
