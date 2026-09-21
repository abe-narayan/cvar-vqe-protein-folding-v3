#!/usr/bin/env python
"""s32/s32_D3_mover.py -- S32 lane D rung D3-M: PHYSICS AS MOVER, not as ranker.

Registered in s32/PREREG_S32_D.md, commit 34973b1b.

WHY THIS AND NOT A SCORE.  Every SCALAR a force field can produce about one structure is a
distance-map reading and therefore adds nothing: energies and components (S31 Cor. B1),
ensemble averages and every dynamical observable (Theorem D-E), normal-mode and curvature
spectra (S31 Cor. B1'), and every contraction of two equivariant fields such as
force-vs-prior-gradient alignment (Theorem D-F).  What is NOT closed is the DISPLACEMENT
itself, used as a correction: `Delta(x) = M(x) - x` is O(3)-EQUIVARIANT, and G1 bounds
invariant scalars only (Corollary D-G).  It is also the shape of the prize: S31 section 20.1
prices a DIRECTION and says it "does not have to be accurate -- it has to point the right
way".

WHAT IS MEASURED, per target, on the production built chain:
  (a) ORACLE / NOT DEPLOYABLE diagnostic -- does the physical move point at the error?
        e_prod   = d(x_prod)      - d_nat      the production chain's own pair-distance error
        e_pool75 = pool75_mean    - d_nat      the pool's common mode (S31's `mu`, recomputed here)
        e_disto  = dg["expected"] - d_nat      the prior's error       (S31's `y`,  recomputed here)
      cos(Delta_d, -e_*) for each.  A helpful move points AGAINST the error, hence the minus.
      Signs and references were fixed in the prereg before any number existed.
  (b) NATIVE-FREE endpoint -- built-chain CA-RMSD of the relaxed chain against the production
      chain, paired, frozen folds, `s24.stats_lib.compare` (LOWER IS BETTER).
  (c) CONTROLS -- (i) a per-target magnitude-matched random direction in the same
      pair-distance space, matched to this arm's OWN ||Delta_d|| (contract rule 6);
      (ii) the G1-closed twin ||Delta_d|| recorded so the directional claim can be priced
      against the closed magnitude; (iii) geometry secondaries on row one (contract rule 15).

BASIS.  Every RMSD here is a BUILT CHAIN.  The production chain is rebuilt in this job from
the shipped top-75 indices, and its mean is asserted against 3.2105 within the measured
per-target projection floor.  No cloud number is differenced against a chain number.

    python s32/s32_D3_mover.py --pass A            # fast restraint ladder
    python s32/s32_D3_mover.py --pass B            # converged, production semantics
    python s32/s32_D3_mover.py --probe 1A13
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I          # noqa: E402
from s24 import stats_lib as ST          # noqa: E402

RES = os.path.join(HERE, "results")

# pass A: the fast restraint ladder.  pass B: production semantics, converged.
LADDER = {"A": [("k100", 100.0, 200, 5.0), ("k10", 10.0, 200, 5.0),
                ("k1", 1.0, 200, 5.0), ("k0", 0.0, 200, 5.0)],
          "B": [("k10c", 10.0, 0, 1.0), ("k1c", 1.0, 0, 1.0)]}


def pair_idx(n, min_sep=2):
    return I.pair_index(n, min_sep=min_sep)


def pdists(X, ii, jj):
    return np.linalg.norm(X[ii] - X[jj], axis=-1)


def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a.dot(b) / (na * nb)) if na > 1e-12 and nb > 1e-12 else np.nan


def virtual_bond(X):
    d = np.linalg.norm(X[1:] - X[:-1], axis=-1)
    return float(d.mean()), float(d.std())


def production(t):
    """Rebuild production for one target: top-75 coordinate average -> STAGE 3b chain."""
    u = I.load_univ(t["pdb"])
    p = I.pool_idx(u, k=500)
    rec = I.shipped_record(t["pdb"])
    sub = np.asarray(rec["sub"], int)
    W = np.asarray(u["W"][p], float)
    C, _ = I.coordinate_average(W[sub])
    pr = I.project(C, t["seq"], int(t["fold"]))
    return {"u": u, "p": p, "sub": sub, "W": W, "C": C,
            "ca": np.asarray(pr["ca"], float),
            "phi": np.asarray(pr["phi"], float), "psi": np.asarray(pr["psi"], float),
            "nat": np.asarray(u["nat_ca"], float)}


def one_target(t, rungs, rng):
    from core import geometry as geo
    from core import amber as ar
    import torsion_lib2 as tl2

    P = production(t)
    n = int(t["n"])
    ii, jj = pair_idx(n)
    nat, ca0 = P["nat"], P["ca"]
    d_nat = pdists(nat, ii, jj)
    d_prod = pdists(ca0, ii, jj)

    # ORACLE reference error vectors, recomputed in this job from this job's own objects.
    e_prod = d_prod - d_nat
    D75 = np.stack([pdists(P["W"][k], ii, jj) for k in P["sub"]])
    e_pool75 = D75.mean(0) - d_nat
    dg = I.distogram(t["pdb"])
    e_disto = np.asarray(dg["expected"], float) - d_nat

    row = {"pdb": t["pdb"], "n": n, "fold": int(t["fold"]),
           "rmsd_prod": float(I.ca_rmsd(ca0, nat)),
           "norm_e_prod": float(np.linalg.norm(e_prod)),
           "norm_e_pool75": float(np.linalg.norm(e_pool75)),
           "norm_e_disto": float(np.linalg.norm(e_disto)),
           "cos_eprod_epool75": cosine(e_prod, e_pool75),
           "cos_eprod_edisto": cosine(e_prod, e_disto)}
    row["vb_mean_prod"], row["vb_sd_prod"] = virtual_bond(ca0)

    BB = geo.build_backbone_batch(P["phi"][None], P["psi"][None])
    cd = {k: v[0] for k, v in BB.items()}
    tab = tl2.library_for(t["seq"], 8, t["seq"])          # target's own sequence held out
    rep = tl2.PerResidueTorsion(t["seq"], tab, chi_bits=False)

    for name, k, steps, tol in rungs:
        t0 = time.time()
        try:
            r = ar.refine_coords(t["seq"], rep, cd, k_restraint=k, steps=steps,
                                 tolerance=tol, components=True, memo=False)
        except Exception as exc:                                   # noqa: BLE001
            row["%s_error" % name] = repr(exc)[:200]
            continue
        finally:
            ar.clear_cache()
        ca1 = np.asarray(r["ca"], float)
        d1 = pdists(ca1, ii, jj)
        dd = d1 - d_prod                                           # the MOVE, in pair space
        row["%s_wall" % name] = round(time.time() - t0, 2)
        row["%s_e0" % name] = float(r["energy_initial"])
        row["%s_e1" % name] = float(r["energy"])
        row["%s_converged" % name] = bool(r["converged"])
        row["%s_rmsd" % name] = float(I.ca_rmsd(ca1, nat))         # BUILT CHAIN, same basis
        row["%s_move_ca" % name] = float(I.ca_rmsd(ca1, ca0))      # rule 15: displacement
        row["%s_norm_dd" % name] = float(np.linalg.norm(dd))       # the G1-CLOSED twin
        row["%s_vb_mean" % name], row["%s_vb_sd" % name] = virtual_bond(ca1)
        row["%s_cos_eprod" % name] = cosine(dd, -e_prod)
        row["%s_cos_epool75" % name] = cosine(dd, -e_pool75)
        row["%s_cos_edisto" % name] = cosine(dd, -e_disto)
        # --- the SCALE ladder.  S31 section 20.1's shrink curve says a correction does not have
        # to be accurate, it has to point the right way, and saturates by c = 0.75.  If the move
        # carries a useful component under a larger useless one, a partial step finds it.  Each
        # alpha is its OWN registered arm; no per-target maximum is taken (order statistic).
        # BASIS: an interpolate of two valid chains, NOT reprojected -- the virtual-bond mean and
        # sd ship on the same row (contract rule 15) so the reader can price the deviation.
        for al in (0.25, 0.50, 0.75):
            xa = ca0 + al * (ca1 - ca0)
            row["%s_a%02d_rmsd" % (name, int(al * 100))] = float(I.ca_rmsd(xa, nat))
            vbm, vbs = virtual_bond(xa)
            row["%s_a%02d_vb_mean" % (name, int(al * 100))] = vbm
            row["%s_a%02d_vb_sd" % (name, int(al * 100))] = vbs
        # CONTROL (i): magnitude-matched random direction in the SAME space, matched to THIS
        # arm's own ||dd||.  Own distribution, not a best draw (contract rule 10): 32 draws.
        cs = []
        for _ in range(32):
            g = rng.standard_normal(len(dd))
            g *= np.linalg.norm(dd) / np.linalg.norm(g)
            cs.append(cosine(g, -e_prod))
        row["%s_ctrl_cos_mean" % name] = float(np.mean(cs))
        row["%s_ctrl_cos_sd" % name] = float(np.std(cs))
        row["%s_ctrl_cos_max" % name] = float(np.max(cs))
    return row


def run(which, shard, of):
    rungs = LADDER[which]
    tg = I.targets()
    mine = [t for i, t in enumerate(tg) if i % of == shard]
    os.makedirs(RES, exist_ok=True)
    jl = os.path.join(RES, "s32_D3_mover_%s_%d_%d.jsonl" % (which, shard, of))
    done = set()
    if os.path.exists(jl):
        with open(jl, encoding="utf-8") as fh:
            for ln in fh:
                try:
                    done.add(json.loads(ln)["pdb"])
                except Exception:                                  # noqa: BLE001
                    pass
    rng = np.random.default_rng(20320032 + shard)
    t0 = time.time()
    for i, t in enumerate(mine):
        if t["pdb"] in done:
            continue
        row = one_target(t, rungs, rng)
        with open(jl, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(row) + "\n")
        print("[%s %d/%d] %s  %.1fs  elapsed %.0fs" % (
            which, i + 1, len(mine), t["pdb"],
            sum(v for k, v in row.items() if k.endswith("_wall") and isinstance(v, float)),
            time.time() - t0), flush=True)
    print("done", jl)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pass", dest="which", default="A", choices=list(LADDER))
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--probe", default=None)
    a = ap.parse_args()
    if a.probe:
        t = {x["pdb"]: x for x in I.targets()}[a.probe]
        r = one_target(t, LADDER[a.which], np.random.default_rng(0))
        print(json.dumps(r, indent=2, default=float))
        return
    run(a.which, a.shard, a.of)


if __name__ == "__main__":
    main()
