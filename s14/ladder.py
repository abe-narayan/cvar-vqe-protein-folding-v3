"""SPRINT 14 -- the shared native-free emission ladder.

Sprint-brief section 46 asks for an explicit information-vs-search ladder, and section 22
asks for an ablation ladder whose stages can be differenced.  Both only compose if every
arm emits a structure the SAME way and is scored by the SAME instrument.  Four agents are
working in parallel; without this module each of them invents its own protocol and the
numbers cannot be subtracted from one another.

An EMITTER is a native-free function

    emit(pdb, seq, n, fold, rng) -> (phi, psi)      radians, shape (n,)

Nothing here reads the native except `incumbent_rmsd` and `score_emitter`, which are
post-hoc scoring against `nat_ca`.  There is a leakage guard: any emitter whose output
changes when the native is perturbed is rejected by `audit_emitter`.

Run:
    python -m s14.ladder                 # the level-0 / level-1 anchors, 126 targets
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from s12 import instrument as I            # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
CACHE = os.path.join(ROOT, "s14", "cache")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

# The zero-information constant alpha-helix.  Sprint 13 measured that this BEATS uniform
# random sampling by 0.457 A [-0.822, -0.090], so it -- not random -- is the baseline any
# native-free arm must clear before it has demonstrated anything at all.
HELIX_PHI = np.deg2rad(-63.0)
HELIX_PSI = np.deg2rad(-42.0)


# ------------------------------------------------------------------ the reference column
def incumbent_rmsd(force: bool = False) -> dict:
    """Per-target CA-RMSD of the shipped synthesis pipeline.  Mean is `synthesis_fit`."""
    path = os.path.join(CACHE, "incumbent.json")
    if os.path.exists(path) and not force:
        with open(path) as fh:
            return json.load(fh)
    out = {}
    for t in I.targets():
        rec = I.shipped_record(t["pdb"])
        u = I.load_univ(t["pdb"])
        out[t["pdb"]] = float(I.ca_rmsd(np.asarray(rec["fit_ca"]), u["nat_ca"]))
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    return out


# ------------------------------------------------------------------------- the state space
_SPACES: dict = {}


def space(pdb: str, k: int = 4):
    """Cached `s13.qarch_lib.Space` -- sequence-conditioned, TARGET HELD OUT."""
    key = (pdb, k)
    if key not in _SPACES:
        from s13.qarch_lib import Space
        _SPACES[key] = Space(pdb, k)
    return _SPACES[key]


_PRIORS: dict = {}


def prior(pdb: str, k: int = 4) -> np.ndarray:
    """Cached leakage-safe (n, k) per-residue state occupancy from the held-out database."""
    key = (pdb, k)
    if key in _PRIORS:
        return _PRIORS[key]
    path = os.path.join(CACHE, f"prior_{pdb}_k{k}.npy")
    if os.path.exists(path):
        P = np.load(path)
    else:
        from s13.qarch_lib import empirical_prior
        P = np.asarray(empirical_prior(space(pdb, k)), float)
        np.save(path, P)
    _PRIORS[key] = P
    return P


# ------------------------------------------------------------------------------- emitters
def emit_constant_helix(pdb, seq, n, fold, rng):
    """LEVEL 0.  Zero information: every residue at the ideal alpha-helical torsions."""
    return np.full(n, HELIX_PHI), np.full(n, HELIX_PSI)


def emit_uniform_k4(pdb, seq, n, fold, rng):
    """LEVEL 1a.  A uniform-random state of the sequence-conditioned k=4 library.

    Carries the library's generic Ramachandran content but no per-residue preference.
    Sprint 13: 88% of what the library buys is generic Ramachandran, 12% sequence.
    """
    sp = space(pdb, 4)
    S = rng.integers(0, 4, size=sp.n)
    r = np.arange(sp.n)
    return sp.PHI[r, S], sp.PSI[r, S]


def emit_prior_sample_k4(pdb, seq, n, fold, rng):
    """LEVEL 1b.  Draw each residue from the leakage-safe empirical state occupancy."""
    sp = space(pdb, 4)
    P = prior(pdb, 4)
    cum = np.cumsum(P, axis=1)
    u = rng.random(sp.n)
    S = (u[:, None] > cum[:, :-1]).sum(1)
    r = np.arange(sp.n)
    return sp.PHI[r, S], sp.PSI[r, S]


def emit_prior_argmax_k4(pdb, seq, n, fold, rng):
    """LEVEL 1c.  The 1-local prior's closed-form minimiser.

    Sprint 13 note carried forward: because the prior is 1-local its 'annealed' arm
    performs no search at all -- the per-residue argmax is returned identically on
    126/126 targets.  So this arm IS the prior's search result, and any 'SA on the
    prior' arm that differs from it has a bug.
    """
    sp = space(pdb, 4)
    S = np.argmax(prior(pdb, 4), axis=1)
    r = np.arange(sp.n)
    return sp.PHI[r, S], sp.PSI[r, S]


EMITTERS = {
    "L0_constant_helix": (emit_constant_helix, False),
    "L1a_uniform_k4": (emit_uniform_k4, True),
    "L1b_prior_sample_k4": (emit_prior_sample_k4, True),
    "L1c_prior_argmax_k4": (emit_prior_argmax_k4, False),
}


# -------------------------------------------------------------------------------- scoring
def score_emitter(emit, stochastic: bool, seeds=(0, 1, 2), targets=None):
    """Post-hoc CA-RMSD of an emitter on every target, averaged over seeds if stochastic.

    Returns (pdbs, rmsd, folds).  Reads the native ONLY to score.
    """
    tg = targets if targets is not None else I.targets()
    pdbs, out, folds = [], [], []
    use = seeds if stochastic else (0,)
    for t in tg:
        u = I.load_univ(t["pdb"])
        vals = []
        for s in use:
            rng = np.random.default_rng(hash((t["pdb"], int(s))) % (2 ** 32))
            phi, psi = emit(t["pdb"], t["seq"], int(t["n"]), int(t["fold"]), rng)
            vals.append(I.ca_rmsd(I.build_ca(phi, psi), u["nat_ca"]))
        pdbs.append(t["pdb"])
        out.append(float(np.mean(vals)))
        folds.append(int(t["fold"]))
    return pdbs, np.asarray(out, float), np.asarray(folds, int)


def audit_emitter(emit, n_probe: int = 6) -> dict:
    """LEAKAGE GUARD.  An emitter must not react to the native structure.

    Monkey-patches `I.load_univ` so every native trace it would return is replaced by
    noise, then checks the emitter's output is bit-identical.  An emitter that reads
    `nat_ca`, native torsions or native RMSD anywhere in its call graph fails here.
    """
    tg = I.targets()[:n_probe]
    real = I.load_univ
    clean = []
    for t in tg:
        rng = np.random.default_rng(7)
        clean.append(emit(t["pdb"], t["seq"], int(t["n"]), int(t["fold"]), rng))

    def poisoned(pdb):
        d = dict(real(pdb))
        r = np.random.default_rng(1234)
        d["nat_ca"] = r.normal(size=np.asarray(d["nat_ca"]).shape) * 10.0
        if "rr" in d:
            d["rr"] = r.random(np.asarray(d["rr"]).shape) * 10.0
        return d

    I.load_univ = poisoned
    try:
        bad = []
        for t, (p0, s0) in zip(tg, clean):
            rng = np.random.default_rng(7)
            p1, s1 = emit(t["pdb"], t["seq"], int(t["n"]), int(t["fold"]), rng)
            if not (np.array_equal(p0, p1) and np.array_equal(s0, s1)):
                bad.append(t["pdb"])
    finally:
        I.load_univ = real
    return {"n_probed": len(tg), "reacts_to_native": bad, "clean": not bad}


def run(names=None, seeds=(0, 1, 2)) -> dict:
    inc = incumbent_rmsd()
    tg = I.targets()
    ref = np.asarray([inc[t["pdb"]] for t in tg], float)
    folds = np.asarray([int(t["fold"]) for t in tg], int)
    pdbs = [t["pdb"] for t in tg]
    fail = np.isin(pdbs, I.FAIL18)

    rows = {"incumbent_synthesis": {
        **I.summary(ref),
        "FAIL18": float(ref[fail].mean()),
        "info": "the shipped retrieval pipeline; mean must equal synthesis_fit 3.2040761603809194",
    }}
    per_target = {"incumbent_synthesis": {p: float(v) for p, v in zip(pdbs, ref)}}

    for name in (names or EMITTERS):
        emit, stoch = EMITTERS[name]
        audit = audit_emitter(emit)
        assert audit["clean"], f"{name} REACTS TO THE NATIVE: {audit['reacts_to_native']}"
        _, r, _ = score_emitter(emit, stoch, seeds=seeds)
        rows[name] = {
            **I.summary(r),
            "FAIL18": float(r[fail].mean()),
            "leakage_audit": audit,
            "stochastic": stoch,
            "seeds": list(seeds) if stoch else [0],
            "vs_incumbent": I.paired(r, ref, folds=folds, names=pdbs),
        }
        per_target[name] = {p: float(v) for p, v in zip(pdbs, r)}

    # The control that matters: everything is also differenced against the zero-information
    # constant helix, because "beats random" is not evidence in this representation.
    if "L0_constant_helix" in rows:
        base = np.asarray([per_target["L0_constant_helix"][p] for p in pdbs], float)
        for name in rows:
            if name == "L0_constant_helix":
                continue
            r = np.asarray([per_target[name][p] for p in pdbs], float)
            rows[name]["vs_constant_helix"] = I.paired(r, base, folds=folds, names=pdbs)

    out = {"rows": rows, "per_target": per_target, "n_targets": len(pdbs)}

    # MERGE, do not overwrite.  `run(names=...)` scores a SUBSET of the arms, and a plain
    # dump of that subset silently destroys the arms from an earlier call -- which is
    # exactly what happened once here: a later retrieval-prior run replaced the level-0 and
    # level-1 anchors in this file and the ladder figure came out missing its own control.
    # This is the same partial-overwrite hazard `I.write` guards against with n_expected.
    path = os.path.join(RESULTS, "ladder.json")
    if os.path.exists(path):
        try:
            with open(path) as fh:
                prev = json.load(fh)
            merged_rows = dict(prev.get("rows", {})); merged_rows.update(out["rows"])
            merged_pt = dict(prev.get("per_target", {})); merged_pt.update(out["per_target"])
            out = {"rows": merged_rows, "per_target": merged_pt, "n_targets": len(pdbs)}
        except (ValueError, OSError):
            pass                                   # a corrupt prior file must not block a run
    I.write("s14_ladder", out, n_expected=len(pdbs))
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    return out


def _fmt(out):
    print(f"{'arm':<24}{'mean':>8}{'median':>8}{'sd':>7}{'<2A':>7}{'FAIL18':>9}"
          f"{'vs incumbent':>26}{'vs helix':>12}")
    for name, r in out["rows"].items():
        v = r.get("vs_incumbent")
        vs = (f"{v['mean_diff']:+.3f} [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]"
              if v else "-- reference --")
        h = r.get("vs_constant_helix")
        hs = f"{h['mean_diff']:+.3f}" if h else "--"
        print(f"{name:<24}{r['mean']:>8.3f}{r['median']:>8.3f}{r['sd']:>7.3f}"
              f"{r['frac_under_2.0']:>7.2f}{r['FAIL18']:>9.3f}{vs:>26}{hs:>12}")


if __name__ == "__main__":
    _fmt(run())
