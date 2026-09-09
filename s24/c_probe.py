"""s24/c_probe.py -- LANE C EXPLORATORY PROBE.  NOT A CLAIM.  NOT PRE-REGISTERED AS DIRECTIONAL.

Purpose: get an honest early read on whether `p(phi,psi | sequence, distogram) -> sample -> build`
can plausibly reach the L2 spec (<= ~3.9 A standalone through the shipped readout), BEFORE any
model is trained and before Lane A's corpus exists.  Run on a small fixed target subset; every
number here is EXPLORATORY and is superseded by the pre-registered n=126 run.

The ladder is ordered by how much information each torsion distribution is allowed to see.
Every arm emits N built chains, is scored by the SAME shipped Bayes-risk distogram functional,
and passes through the SAME readout (top-75, uniform coordinate average in the medoid frame).
Members are BUILT CHAINS; the emitted average is a POINT CLOUD (an average of 75 chains is not
itself an ideal-geometry chain), so it is commensurable with the incumbent's 3.0483 A point cloud.

  R  REBUILD-75   the incumbent's own top-75 windows, rebuilt from THEIR OWN torsions.
                  Not a generator -- the calibration that prices the ideal-geometry manifold.
  T0 HELIX        constant alpha-helix (-63, -42) + jitter.  The plausible zero-information
                  control (uniform-on-the-torus is a worse measure, not an uninformative one).
  T1 BLIND        sequence-blind Ramachandran marginal of this target's own legal universe.
  T2 RESTYPE      per-residue-type Ramachandran marginal of the same universe.  This is the
                  measured ceiling of the sequence->torsion channel (36.9/68.9 deg).
  T3 POOL         per-residue 2-component von Mises mixture fit to the top-75 RETRIEVED windows'
                  torsions (s19 qb_lib.fit_basins).  s14 measured this channel at phi 33.6 /
                  psi 59.2 deg -- BETTER than the trained leave-fold-out sequence predictor's
                  36.1 / 62.4.  It is therefore an ACHIEVABLE UPPER REFERENCE for any learned
                  p(phi,psi | sequence, distogram): a network has to beat a free construction
                  that already beats the best trained torsion predictor this project has built.

ORACLE labelling: every CA-RMSD reads the native and is post-hoc scoring of a native-free
decision.  Nothing here selects on the native.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

from s12 import instrument as I          # noqa: E402
from s15 import seed as SD               # noqa: E402
from s19 import qb_lib as QB             # noqa: E402
from core import project as pj           # noqa: E402

TOPM = 75
NSAMP = 2000
HELIX = (np.deg2rad(-63.0), np.deg2rad(-42.0))


def wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


# --------------------------------------------------------------------------- samplers
def s_helix(n, N, rng, jitter=np.deg2rad(15.0)):
    phi = HELIX[0] + jitter * rng.standard_normal((N, n))
    psi = HELIX[1] + jitter * rng.standard_normal((N, n))
    return wrap(phi), wrap(psi)


def s_blind(u, n, N, rng):
    """Sequence-blind: resample (phi,psi) PAIRS from the whole legal universe, i.i.d. per residue.

    Pairs, not independent marginals -- an independent phi/psi marginal is not a Ramachandran
    distribution and would be a strawman.
    """
    P = wrap(u["PHI"]).ravel(); S = wrap(u["PSI"]).ravel()
    k = rng.integers(0, len(P), (N, n))
    return P[k], S[k]


def s_restype(u, n, N, rng):
    """Per-residue-type Ramachandran marginal of the same legal universe."""
    P = wrap(u["PHI"]).ravel(); S = wrap(u["PSI"]).ravel()
    C = np.asarray(u["S"], int).ravel()
    by = {}
    for c in range(20):
        w = np.where(C == c)[0]
        by[c] = w
    tgt_codes = [I.ALPHABET.index(a) for a in u["seq"]]
    phi = np.zeros((N, n)); psi = np.zeros((N, n))
    for i in range(n):
        pool = by[tgt_codes[i]]
        if len(pool) < 50:                      # fall back to the blind marginal
            pool = np.arange(len(P))
        k = pool[rng.integers(0, len(pool), N)]
        phi[:, i] = P[k]; psi[:, i] = S[k]
    return phi, psi


def s_pool(u, n, N, rng, idx_top):
    """s19 fit_basins on the top-75 RETRIEVED windows' torsions, then draw."""
    PHI = wrap(u["PHI"][idx_top]); PSI = wrap(u["PSI"][idx_top])
    mu, kap, w = QB.fit_basins(PHI, PSI, rng)
    bits = (rng.random((N, n)) > w[None, :, 0]).astype(np.int64)
    return QB.draw_from_basins(bits, mu, kap, rng)


# --------------------------------------------------------------------------- readout
def _avg(members):
    P = I.pairwise_rmsd(members); b = I.medoid(P)
    return I.superpose_batch(members, members[b]).mean(0)


def _kabsch_R(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def _bias(C, nat):
    R = _kabsch_R(C, nat)
    return (C - C.mean(0)) @ R.T - (nat - nat.mean(0))


def _cos(a, b):
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    return float((a * b).sum() / (na * nb)) if na > 0 and nb > 0 else float("nan")


def readout(CA, dg, ij, nat, m=TOPM, rng=None):
    """Two readouts of the SAME emitted set, to separate SOURCE bias from SCORE bias.

    `scored`  top-m by the shipped Bayes-risk functional -- the incumbent's own selector, so
              any bias the SELECTOR imposes is shared with the incumbent by construction.
    `uniform` m taken uniformly from the emitted set -- no scoring at all.  This is the readout
              L2's library arm actually used (its cos 0.647 is an UNSCORED number), so it is the
              only one commensurable with the published 0.65 bar.
    """
    i, j = ij
    sc = np.asarray(I.shipped_score(dg, I.pair_dists(CA, i, j)), float)
    o = np.argsort(sc, kind="stable")[:m]
    sel = CA[o]
    c = _avg(sel)
    ru = (rng or np.random.default_rng(0)).choice(len(CA), m, replace=False)
    cu = _avg(CA[ru])
    rr = I.kabsch_rmsd_batch(CA, nat)                       # ORACLE, post-hoc
    return {"rmsd": float(I.ca_rmsd(c, nat)), "cloud": c,
            "rmsd_unif": float(I.ca_rmsd(cu, nat)), "cloud_unif": cu,
            "member_mean_ORACLE": float(I.kabsch_rmsd_batch(sel, nat).mean()),
            "gen_best_ORACLE": float(rr.min()), "gen_mean_ORACLE": float(rr.mean())}


def run(pdbs=None, nsamp=NSAMP):
    tg = I.targets()
    if pdbs is None:
        pdbs = [t["pdb"] for t in tg[::16]]                  # a fixed, arbitrary-by-index subset
    sel = [t for t in tg if t["pdb"] in set(pdbs)]
    rows = []
    for t in sel:
        pdb = t["pdb"]; n = t["n"]
        u = I.load_univ(pdb); nat = u["nat_ca"]
        dg = I.distogram(pdb, u["seq"], u["fold"]); ij = I.pair_index(n)
        idx = I.pool_idx(u); Wp = u["W"][idx]
        scp = np.asarray(I.shipped_score(dg, I.pair_dists(Wp, *ij)), float)
        top = np.argsort(scp, kind="stable")[:TOPM]

        # --- incumbent (real window coordinates, point cloud)
        cA = _avg(Wp[top]); inc = float(I.ca_rmsd(cA, nat)); eA = _bias(cA, nat)

        # --- R: rebuild the SAME 75 windows from their own torsions
        CAr = np.asarray(pj.build_ca_exact(wrap(u["PHI"][idx][top]), wrap(u["PSI"][idx][top])), float)
        cR = _avg(CAr); reb = float(I.ca_rmsd(cR, nat))

        r = {"pdb": pdb, "n": n, "fold": int(u["fold"]), "incumbent": inc, "rebuild75": reb,
             "cos_rebuild": _cos(_bias(cR, nat), eA)}
        rng0 = SD.stable_rng("c_probe", pdb)
        for tag, fn in (("T0_helix", lambda g: s_helix(n, nsamp, g)),
                        ("T1_blind", lambda g: s_blind(u, n, nsamp, g)),
                        ("T2_restype", lambda g: s_restype(u, n, nsamp, g)),
                        ("T3_pool", lambda g: s_pool(u, n, nsamp, g, idx[top]))):
            ph, ps = fn(rng0)
            CA = np.asarray(pj.build_ca_exact(ph, ps), float)
            o = readout(CA, dg, ij, nat, rng=rng0)
            r[tag] = {k: v for k, v in o.items() if not k.startswith("cloud")}
            r[tag]["cos_vs_incumbent"] = _cos(_bias(o["cloud"], nat), eA)
            r[tag]["cos_unif_vs_incumbent"] = _cos(_bias(o["cloud_unif"], nat), eA)
            # within-source control: an independent draw from the SAME sampler
            ph2, ps2 = fn(rng0)
            CA2 = np.asarray(pj.build_ca_exact(ph2, ps2), float)
            o2 = readout(CA2, dg, ij, nat, rng=rng0)
            r[tag]["cos_within"] = _cos(_bias(o["cloud"], nat), _bias(o2["cloud"], nat))
            r[tag]["cos_unif_within"] = _cos(_bias(o["cloud_unif"], nat), _bias(o2["cloud_unif"], nat))
        rows.append(r)
        print("  %s n=%d inc %.3f reb %.3f | %s" % (
            pdb, n, inc, reb,
            "  ".join("%s %.3f/c%.2f" % (k[:2], r[k]["rmsd"], r[k]["cos_vs_incumbent"])
                      for k in ("T0_helix", "T1_blind", "T2_restype", "T3_pool"))), flush=True)
    p = os.path.join(RES, "c_probe.json"); tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"rows": rows, "nsamp": nsamp, "EXPLORATORY": True}, fh)
    os.replace(tmp, p)
    rep(rows)
    return rows


def rep(rows):
    g = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    print("\nEXPLORATORY, n=%d targets, N=%d samples/target, top-%d uniform average, POINT CLOUD."
          % (len(rows), NSAMP, TOPM))
    print("            |------ SCORED top-75 (the incumbent's OWN selector) ------|"
          "  |---- UNIFORM 75 (no score) ----|")
    print("  %-12s %8s %8s %8s %7s %8s | %8s %8s %8s %7s"
          % ("arm", "RMSD", "cos_inc", "cos_wtn", "ratio", "genbest",
             "RMSD", "cos_inc", "cos_wtn", "ratio"))
    print("  %-12s %8.4f" % ("INCUMBENT", g("incumbent").mean()))
    print("  %-12s %8.4f %8.4f" % ("REBUILD75", g("rebuild75").mean(), g("cos_rebuild").mean()))
    for k in ("T0_helix", "T1_blind", "T2_restype", "T3_pool"):
        f = lambda q: np.array([r[k][q] for r in rows], float)      # noqa: E731
        print("  %-12s %8.4f %8.4f %8.4f %7.3f %8.3f | %8.4f %8.4f %8.4f %7.3f"
              % (k, f("rmsd").mean(), f("cos_vs_incumbent").mean(), f("cos_within").mean(),
                 f("cos_vs_incumbent").mean() / f("cos_within").mean(), f("gen_best_ORACLE").mean(),
                 f("rmsd_unif").mean(), f("cos_unif_vs_incumbent").mean(),
                 f("cos_unif_within").mean(),
                 f("cos_unif_vs_incumbent").mean() / f("cos_unif_within").mean()))
    print("  SPEC: standalone <= ~3.9 A AND cos <= ~0.65 merely to break even.")
    print("  L2's published 0.647 / 0.933 / ratio 0.693 was measured on an UNSCORED library draw.")


if __name__ == "__main__":
    run()
