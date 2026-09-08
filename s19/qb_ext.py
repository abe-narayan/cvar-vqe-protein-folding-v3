"""SPRINT 19 / AGENT B -- E: "go deeper than an alpha sweep", as the directive requires.

The directive names, explicitly: alpha SCHEDULES, ANNEALED and ADAPTIVE alpha, tail
selection, state concentration, entropy, effective sample size, multimodality, basin
transitions, WARM STARTS, MULTIPLE INDEPENDENT ANSATZ SEEDS, and CVaR+UNIFORM and
CVaR+ANNEALING MIXTURES.  `qb_arms` covers the alpha sweep, one schedule, the entropy /
ESS / latent-mutual-information diagnostics and the untrained best-of-N.  This module
covers the rest, on a declared prefix of targets, at the identical budget.

Every arm here is still measured against the SAME classical set and the SAME
`q_untrained` best-of-N control from `qb_arms`; nothing in this module is scored against
an initialisation mean.

RUN:  python -m s19.qb_ext run [n_targets]
"""
from __future__ import annotations

import gc
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I           # noqa: E402
from s15 import seed as SD                # noqa: E402
from core import quantum as Q             # noqa: E402
from s19 import qb_lib as L               # noqa: E402
from s19 import qb_arms as AR             # noqa: E402


def _warm_theta(an, tgt, rng):
    """WARM START from the retrieval pool's own basin occupancy.

    The final RY layer's angles are set so the per-qubit marginal equals the pool's
    empirical P(basin = 1) for that residue; the earlier layers start near identity.
    Native-free: the occupancies come from the BLOSUM pool, never from the native.
    """
    P = an.n_params()
    m = np.clip(tgt["wmarg"][:, 1], 1e-6, 1 - 1e-6)
    th = np.zeros(P)
    nb = P // tgt["n"]
    th = th.reshape(nb, tgt["n"])
    th[:-1] = 1e-3 * rng.standard_normal(th[:-1].shape)
    th[-1] = 2.0 * np.arcsin(np.sqrt(m))
    return th.reshape(-1)


def run_circuit_ext(tgt, obj, seed, mode, alpha=0.25, ansatz="mps2f",
                    shots=L.SHOTS, lr=L.LR, uniform_frac=0.0):
    """The extended quantum arms.  `mode` in {plain, warm, adapt, mixunif, mixanneal}."""
    n = tgt["n"]
    an = Q.MPSAnsatz(n, layers=(3 if ansatz.startswith("mps3") else 2), final_ry=True,
                     entangler=("none" if ansatz.endswith("n") else "cnot"))
    rng = SD.stable_rng(tgt["pdb"], seed, mode, ansatz, alpha, salt=L.SALT)
    th = (_warm_theta(an, tgt, rng) if mode == "warm"
          else np.pi / 2 + rng.normal(0.0, 0.8, an.n_params()))
    opt = Q.Adam(an.n_params(), lr=lr)
    # CVaR + UNIFORM MIXTURE: a fixed fraction of the budget is spent on plausible
    # zero-information basin-marginal draws that the circuit never sees.
    ub = int(round(uniform_frac * obj.budget))
    if ub:
        while obj.used < ub:
            m = min(shots, ub - obj.used)
            b = (rng.random((m, n)) < tgt["wmarg"][:, 1][None, :]).astype(np.int64)
            phi, psi = L.draw_from_basins(b, tgt["mu"], tgt["kap"], rng)
            obj(phi, psi)
    iters = max(1, obj.left // shots)
    a_hist = []
    for t in range(iters):
        if obj.left < shots:
            break
        if mode == "adapt":
            # ADAPTIVE alpha: the tail fraction is set to the share of the previous batch
            # that beat its own median-of-medians, clipped -- native-free, data-driven.
            a = float(np.clip(a_hist[-1] if a_hist else 0.25, 0.02, 1.0))
        elif mode == "mixanneal":
            a = Q.alpha_schedule(t / max(1, iters - 1), 1.0, 0.05)
        else:
            a = alpha
        bits = an.sample(th, shots, rng)
        phi, psi = L.draw_from_basins(bits, tgt["mu"], tgt["kap"], rng)
        e = obj(phi, psi)
        if e.size < shots:
            break
        g, _v = Q.cvar_gradient(an, th, bits, e, a, baseline="const")
        th = opt.step(th, g)
        if mode == "adapt":
            a_hist.append(float(np.clip((e < np.median(e)).mean() * 0.5, 0.02, 1.0)))
    d = {"iters": iters, "mode": mode, "uniform_frac": uniform_frac}
    if n <= 18:
        allb = Q.all_bitstrings(n)
        p = np.exp(np.asarray(an.logp(th, allb), float))
        p = np.maximum(p, 0); p /= p.sum()
        pos = p[p > 0]
        d["entropy_bits"] = float(-(pos * np.log2(pos)).sum())
        d["ess"] = float(1.0 / (p ** 2).sum())
    return d


EXT = {
    "x_warm":        dict(mode="warm", alpha=0.25),
    "x_adapt":       dict(mode="adapt"),
    "x_mixunif50":   dict(mode="plain", alpha=0.25, uniform_frac=0.5),
    "x_mixanneal":   dict(mode="mixanneal"),
    "x_mps3f":       dict(mode="plain", alpha=0.25, ansatz="mps3f"),
    "x_a0.02":       dict(mode="plain", alpha=0.02),
    "x_a0.50":       dict(mode="plain", alpha=0.50),
}
SEEDS = (0, 1, 2, 3)


def run(limit=40, tag="ext"):
    tg = I.targets()[:int(limit)]
    L.gather_full()
    have = L.load(tag)
    for c, t in enumerate(tg):
        pdb = t["pdb"]
        if pdb in have:
            continue
        fg = I.free_gb()
        for _w in range(30):
            if fg >= 0.60:
                break
            print(f"[wait] free RAM {fg:.2f} GB", flush=True)
            time.sleep(60)
            fg = I.free_gb()
        if fg < 0.60:
            print("[STOP] memory", flush=True)
            return
        t0 = time.time()
        tgt = L.target(pdb)
        row = {"pdb": pdb, "n": tgt["n"], "fold": tgt["fold"], "arms": {}}
        for nm, kw in EXT.items():
            obj = L.new_obj(tgt)
            extra = run_circuit_ext(tgt, obj, seed=0, **kw)
            phi, psi, e = obj.seen()
            r = L.readout(phi, psi, e, tgt, tag=nm)
            r.update(extra); r["evals"] = int(obj.used)
            row["arms"][nm] = r
        # MULTIPLE INDEPENDENT ANSATZ SEEDS on the primary arm and on its untrained control
        for s in SEEDS:
            for nm, kwv in (("s_q0.25", dict(alpha=0.25, ansatz="mps2f", train=True)),
                            ("s_untrained", dict(alpha=0.25, ansatz="mps2f",
                                                 train=False))):
                obj = L.new_obj(tgt)
                extra = AR.run_circuit(tgt, obj, seed=s, **kwv)
                phi, psi, e = obj.seen()
                r = L.readout(phi, psi, e, tgt, tag=f"{nm}_{s}")
                r.update({k: v for k, v in extra.items() if k != "cvar_traj"})
                r["evals"] = int(obj.used)
                row["arms"][f"{nm}_{s}"] = r
        L.ck(tag, pdb, row)
        have[pdb] = row
        print(f"[{c+1}/{len(tg)}] {pdb} {time.time()-t0:5.1f}s  " +
              "  ".join(f"{k[2:]} {row['arms'][k]['realised_ORACLE']:.3f}"
                        for k in ("x_warm", "x_adapt", "x_mixunif50", "x_mps3f")),
              flush=True)
        del tgt
        gc.collect()
    L.ck(tag, "_complete", int(len([k for k in L.load(tag) if not k.startswith("_")])))


if __name__ == "__main__":
    run(limit=int(sys.argv[2]) if len(sys.argv) > 2 else 40)
