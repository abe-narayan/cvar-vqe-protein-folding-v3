"""SPRINT 14, coordinator -- the direct chain, end to end, instead of a chained correlation.

WHAT C22 LEFT UNMEASURED, AND WHY IT MATTERS.

C22 measured `native-free proxy -> native Rg` on 126 targets (best r = 0.365,
length-residualised, CI [+0.155, +0.577]) and quoted the OBJ workstream's
`per-target skill <-> native Rg` link at +0.909 from a smaller enumerated set. Composing the
two into "about 61% sign accuracy" assumes linearity and independent errors, and I flagged it
in C22 as an estimate rather than a measurement.

Both halves are now on disk. This measures the chain **directly**:

    does a NATIVE-FREE compactness proxy predict the PER-TARGET SIGN of the learned
    objective's cross-target skill?

That is the actual question the next sprint would ask, and it is answerable today on the 12
targets where both quantities exist.

THREE ARMS, in increasing honesty:

  A  ORACLE ceiling.  native Rg (z-scored) -> per-target cross-target skill.
     This re-derives the +0.909 flip diagnostic in independent code. If it does not reproduce,
     the whole premise is unsafe and that is the finding.

  B  THE CHAIN.  native-free proxy -> per-target cross-target skill, directly.
     No composition, no linearity assumption. This is the number C22 could not supply.

  C  SIGN ACCURACY.  the binary question a conditioning signal actually has to answer:
     does the proxy get the SIGN of the per-target skill right, against a 50% baseline and
     against the majority-class rate (which a constant predictor achieves for free)?

n = 12 is small and the result will be reported as underpowered whatever it says. It is run
because a direct measurement at n=12 is worth more than a composed estimate at n=126, and
because the next sprint should start from a measured number rather than my arithmetic.

ORACLE only in arm A and in the skill labels, which are post-hoc evaluation.

Run:
    python -m s14.signchain
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
from s14 import signpred as SP             # noqa: E402

RESULTS = os.path.join(ROOT, "s14", "results")
CEILING = os.path.join(ROOT, "s12", "results", "s14_obj_ceiling.json")


def load_ceiling():
    with open(CEILING) as fh:
        d = json.load(fh)
    return {r["pdb"]: r for r in d["rows"]}


def load_signpred():
    p = os.path.join(RESULTS, "signpred.json")
    with open(p) as fh:
        d = json.load(fh)
    return {r["pdb"]: r for r in d["per_target"]}


def residualise(x, n_):
    A = np.column_stack([np.ones_like(n_), n_])
    beta, *_ = np.linalg.lstsq(A, x, rcond=None)
    return x - A @ beta


def perm_p(x, y, fn=SP.pearson, n_perm=20000, seed=0):
    """Permutation p-value for |corr|, which is the right test at n=12."""
    rng = np.random.default_rng(seed)
    obs = abs(fn(x, y))
    y = np.asarray(y, float)
    cnt = sum(abs(fn(x, rng.permutation(y))) >= obs for _ in range(n_perm))
    return float((cnt + 1) / (n_perm + 1))


def run(skill_key="cross_rho"):
    ceil = load_ceiling()
    sp = load_signpred()
    pdbs = sorted(set(ceil) & set(sp))

    n = np.asarray([sp[p]["n"] for p in pdbs], float)
    skill = np.asarray([ceil[p][skill_key] for p in pdbs], float)
    rg_nat = np.asarray([sp[p]["rg_native_ORACLE"] for p in pdbs], float)
    proxies = {
        "distogram-predicted Rg": np.asarray([sp[p]["rg_disto"] for p in pdbs], float),
        "retrieval pool mean Rg": np.asarray([sp[p]["rg_pool"] for p in pdbs], float),
        "incumbent emitted Rg": np.asarray([sp[p]["rg_emit"] for p in pdbs], float),
    }

    out = {"n_targets": len(pdbs), "targets": pdbs, "skill_key": skill_key,
           "skill": {p: float(ceil[p][skill_key]) for p in pdbs}, "arms": {}}

    # ---- arm A: ORACLE ceiling -- does the flip diagnostic reproduce in independent code?
    zr = residualise(rg_nat, n)
    out["arms"]["A_ORACLE_native_Rg"] = {
        "pearson": SP.pearson(zr, skill), "spearman": SP.spearman(zr, skill),
        "ci95": SP.boot_ci(zr, skill), "perm_p": perm_p(zr, skill)}

    # ---- arm B: the chain, measured directly
    for name, v in proxies.items():
        vr = residualise(v, n)
        out["arms"]["B_" + name] = {
            "pearson": SP.pearson(vr, skill), "spearman": SP.spearman(vr, skill),
            "ci95": SP.boot_ci(vr, skill), "perm_p": perm_p(vr, skill)}

    # ---- arm C: the binary question a conditioning signal must answer
    s_sign = skill > np.median(skill)
    major = max(s_sign.mean(), 1 - s_sign.mean())
    out["majority_class_rate"] = float(major)
    for name, v in proxies.items():
        vr = residualise(v, n)
        p_sign = vr > np.median(vr)
        acc = float(max((p_sign == s_sign).mean(), (p_sign != s_sign).mean()))
        out["arms"]["C_sign_" + name] = {"accuracy": acc, "vs_majority": acc - major}
    zs = zr > np.median(zr)
    out["arms"]["C_sign_ORACLE_native_Rg"] = {
        "accuracy": float(max((zs == s_sign).mean(), (zs != s_sign).mean())),
        "vs_majority": float(max((zs == s_sign).mean(),
                                 (zs != s_sign).mean())) - major}

    with open(os.path.join(RESULTS, "signchain.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    I.write("s14_signchain", out, n_expected=len(pdbs))

    print(f"THE DIRECT CHAIN: does a compactness proxy predict per-target cross-target "
          f"skill?\nn = {len(pdbs)} targets, skill = {skill_key}\n")
    print(f"{'arm':<42}{'pearson':>9}{'spearman':>10}{'perm p':>9}{'CI95':>22}")
    for k, v in out["arms"].items():
        if not k.startswith("C_"):
            print(f"{k:<42}{v['pearson']:>9.3f}{v['spearman']:>10.3f}{v['perm_p']:>9.3f}"
                  f"  [{v['ci95'][0]:+.3f},{v['ci95'][1]:+.3f}]")
    print(f"\nSIGN ACCURACY (the binary question), majority-class baseline "
          f"{out['majority_class_rate']:.3f}:")
    for k, v in out["arms"].items():
        if k.startswith("C_"):
            print(f"  {k[7:]:<40}{v['accuracy']:>8.3f}   vs majority {v['vs_majority']:+.3f}")
    print("\nReference: the OBJ workstream reported per-target skill correlating +0.909 with "
          "the native Rg z-score.")
    return out


if __name__ == "__main__":
    run()
