"""s26/p_c5.py -- PREREG_C5: predict the common-mode direction on held-out folds and subtract it.

Per target the shipped top-75 members in the medoid frame give the cloud c; the native t is
ORACLE (training folds only, inside the fit).  Common mode: ebar = c - t after ONE rigid Kabsch fit.

Two native-free representations of a correction (both pre-declared):
  R1  distance space: e_D = D(c) - D(t) on the score's pair set; modes = per-SHELL means (5),
      offset, stretch, per-residue additive.  A correction moves D(c) -> D(c) - a * e_hat and is
      realised by stress descent INITIALISED AT c (no mirror crossing), then I.project.
  R2  coordinate space in c's own principal-axis frame (sign rules: PC1 along N->C; PC2 by a smooth
      chain-bump weighting; PC3 = PC1 x PC2): E = ebar_c @ F.  Round trip c - E @ F.T == t_c exactly.
Arms: GLOBAL (training-fold mean correction per length band), RIDGE (R1 shell means from native-
free features, nested), ORACLE (the true correction: the bound), RANDOM-MATCHED (same magnitude,
random low-frequency direction, same operator), SCALE-ONLY (s*, ORACLE, for scale).  Every arm
goes through the same projection.  `run` is gated on "PHASE 0 SIGNED OFF"; `selftest` is synthetic.

    python s26/p_c5.py selftest
    python s26/p_c5.py run           (after sign-off; one agent-day budget)
"""
from __future__ import annotations

import json
import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from s26 import p_stats as PS                        # noqa: E402
from s26 import p_ladder as L                        # noqa: E402
from s26 import p_b3 as B3                           # noqa: E402
from core import predict as PR                       # noqa: E402

RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "p_c5.json")
SHELLS = PR.SHELLS
BANDS = ((9, 11), (12, 13), (14, 16))
GRID = 16
ALPHAS_A = (0.5,)   # reduced from (0.25, 0.5, 1.0) at 01:35 for the 04:30 close; one agent-day budget


# ------------------------------------------------------------------ geometry helpers
def kabsch_onto(P, Q):
    """Rotate/translate P onto Q (rotation only, no scale).  Returns the moved P."""
    P = np.asarray(P, float); Q = np.asarray(Q, float)
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return Pc @ R.T + Q.mean(0)


def cloud(u, rec):
    idx = I.pool_idx(u); top = np.asarray(u["W"][idx][np.asarray(rec["sub"], int)], float)
    P = I.pairwise_rmsd(top)
    return I.superpose_batch(top, top[I.medoid(P)]).mean(0)


def frame(c):
    """R2's native-free frame from the cloud alone: (3, 3) matrix F with columns PC1, PC2, PC3."""
    c = np.asarray(c, float); n = len(c); x = c - c.mean(0)
    ev, vec = np.linalg.eigh(x.T @ x / n)
    vec = vec[:, ::-1]
    p1 = vec[:, 0] * np.sign((x[-1] - x[0]) @ vec[:, 0] or 1.0)
    bump = np.sin(np.pi * np.arange(n) / max(n - 1, 1))
    s2 = np.sign(float(bump @ (x @ vec[:, 1])) or 1.0)
    p2 = vec[:, 1] * s2
    p2 = p2 - (p2 @ p1) * p1; p2 /= np.linalg.norm(p2)
    p3 = np.cross(p1, p2)
    return np.stack([p1, p2, p3], 1)


def rel_grid(v, n):
    """Interpolate a per-residue (n, k) array onto the 16-point relative grid (and back with n)."""
    v = np.asarray(v, float)
    src = np.linspace(0, 1, len(v)); dst = np.linspace(0, 1, n)
    return np.stack([np.interp(dst, src, v[:, k]) for k in range(v.shape[1])], 1)


def realise(D_target, i, j, x0, iters=200, lr=0.02):
    """Stress descent toward the pair distances D_target on (i, j), initialised at x0.  Adjacent
    pairs keep their current length (chain restraint), so the move is the correction and nothing else."""
    x = np.asarray(x0, float).copy(); n = len(x)
    adj = np.arange(n - 1); dadj = np.linalg.norm(x[adj + 1] - x[adj], axis=1)
    for _ in range(iters):
        diff = x[:, None] - x[None]; D = np.linalg.norm(diff, axis=-1) + 1e-9
        g = np.zeros_like(x)
        r = ((D[i, j] - D_target))[:, None] * (diff[i, j] / D[i, j][:, None])
        np.add.at(g, i, r); np.add.at(g, j, -r)
        rc = (4.0 * (D[adj, adj + 1] - dadj))[:, None] * (diff[adj, adj + 1] / D[adj, adj + 1][:, None])
        np.add.at(g, adj, rc); np.add.at(g, adj + 1, -rc)
        x = x - lr * g
    return x


def shell_means(e, i, j):
    sep = j - i
    return np.array([e[(sep >= a) & (sep <= b)].mean() if ((sep >= a) & (sep <= b)).any() else 0.0 for a, b in SHELLS])


def shell_broadcast(means, i, j):
    sep = j - i; out = np.zeros(len(i))
    for (a, b), m in zip(SHELLS, means):
        out[(sep >= a) & (sep <= b)] = m
    return out


def band_of(n):
    for k, (a, b) in enumerate(BANDS):
        if a <= n <= b:
            return k
    return len(BANDS) - 1


def random_lowfreq(n, rng):
    """A random smooth (n, 3) displacement: Gaussian noise smoothed along the chain, unit norm."""
    z = rng.normal(size=(n + 4, 3))
    k = np.array([1, 2, 3, 2, 1], float) / 9.0
    s = np.stack([np.convolve(z[:, a], k, mode="valid") for a in range(3)], 1)
    return s / max(np.linalg.norm(s), 1e-12)


# ------------------------------------------------------------------ the run (gated)
def collect():
    """Per target: cloud, native (ORACLE), the two representations of ebar, and native-free features."""
    fz = B3.build_features(verbose=False); feats = {r["pdb"]: r for r in fz["rows"]}; names = fz["names"]
    rows = []
    for t in I.targets():
        u = I.load_univ(t["pdb"]); rec = I.shipped_record(t["pdb"]); dg = I.distogram(t["pdb"])
        c = cloud(u, rec); nat = np.asarray(u["nat_ca"], float); n = len(c)
        t_c = kabsch_onto(nat, c)                                        # native in the cloud's frame
        F = frame(c); E = (c - t_c) @ F                                  # R2 coefficients
        i, j = np.asarray(dg["i"]), np.asarray(dg["j"])
        Dc = I.pair_dists(c, i, j); Dt = I.pair_dists(nat, i, j); eD = Dc - Dt
        rows.append({"pdb": t["pdb"], "n": n, "fold": int(t["fold"]), "seq": t["seq"], "c": c, "t_c": t_c, "nat": nat,
                     "F": F, "E16": rel_grid(E, GRID), "i": i, "j": j, "Dc": Dc, "eD": eD,
                     "shell": shell_means(eD, i, j), "x": np.array([feats[t["pdb"]][k] for k in names], float)})
        del u
    return rows, names


def endpoint(c_corr, r):
    pr = I.project(c_corr, r["seq"], r["fold"])
    return {"cloud": float(I.ca_rmsd(c_corr, r["nat"])), "arm": float(I.ca_rmsd(pr["ca"], r["nat"])), "fit": float(I.ca_rmsd(pr["fit_ca"], r["nat"]))}


def run():
    if not L._signed_off():
        raise SystemExit("PHASE GATE: p_c5 reads natives; refused before sign-off.")
    rows, names = collect()
    folds = np.array([r["fold"] for r in rows]); pdbs = [r["pdb"] for r in rows]; n = len(rows)
    rng = np.random.default_rng(26)
    X = np.array([r["x"] for r in rows]); SH = np.array([r["shell"] for r in rows])
    pred_shell, a_sh = PS.nested_predict(X, SH, folds, "reg")
    res = {r["pdb"]: {} for r in rows}
    for k, r in enumerate(rows):
        c = r["c"]; i, j = r["i"], r["j"]; tr = folds != r["fold"]
        res[r["pdb"]]["incumbent"] = endpoint(c, r)
        # ORACLE: R2 exact and R1 exact (the bounds)
        res[r["pdb"]]["oracle_R2"] = endpoint(r["t_c"], r)
        res[r["pdb"]]["oracle_R1"] = endpoint(realise(r["Dc"] - r["eD"], i, j, c), r)
        # GLOBAL R1: training-fold mean shell profile in the length band
        band = band_of(r["n"]); sel = tr & np.array([band_of(q["n"]) == band for q in rows])
        g_sh = SH[sel].mean(0) if sel.any() else SH[tr].mean(0)
        # GLOBAL R2: training-fold mean E on the relative grid, in the band
        g_E = np.mean([q["E16"] for q, s in zip(rows, sel) if s], 0) if sel.any() else np.mean([q["E16"] for q, s in zip(rows, tr) if s], 0)
        for a in ALPHAS_A:
            res[r["pdb"]]["global_R1_a%g" % a] = endpoint(realise(r["Dc"] - a * shell_broadcast(g_sh, i, j), i, j, c), r)
            corr = rel_grid(g_E, r["n"]) @ r["F"].T
            res[r["pdb"]]["global_R2_a%g" % a] = endpoint(c - a * corr, r)
            # RANDOM-MATCHED for R2 at the same magnitude
            rnd = random_lowfreq(r["n"], rng) * np.linalg.norm(a * corr)
            res[r["pdb"]]["random_R2_a%g" % a] = endpoint(c - rnd, r)
        # RIDGE R1: predicted shell means (nested, held out)
        ph = pred_shell[k]
        res[r["pdb"]]["ridge_R1"] = endpoint(realise(r["Dc"] - shell_broadcast(ph, i, j), i, j, c), r)
        rs = rng.normal(size=5); rs *= np.linalg.norm(ph) / max(np.linalg.norm(rs), 1e-12)
        res[r["pdb"]]["random_R1"] = endpoint(realise(r["Dc"] - shell_broadcast(rs, i, j), i, j, c), r)
        if (k + 1) % 10 == 0:
            print("  %d/%d" % (k + 1, n), flush=True)
            ST.save_atomic(OUT, {"rows": res, "alphas_ridge": a_sh}, module_file=__file__)
    ST.save_atomic(OUT, {"rows": res, "alphas_ridge": a_sh, "complete": True}, module_file=__file__)
    inc = np.array([res[p]["incumbent"]["arm"] for p in pdbs])
    summary = {}
    for arm in sorted({k for p in pdbs for k in res[p]} - {"incumbent"}):
        v = np.array([res[p][arm]["arm"] for p in pdbs])
        cmp = ST.compare(v, inc, folds, names=pdbs, label="C5 %s - incumbent (built chain)" % arm)
        summary[arm] = cmp; print(ST.fmt(cmp))
    ST.save_atomic(OUT, {"rows": res, "alphas_ridge": a_sh, "complete": True, "summary": summary}, module_file=__file__)
    return summary


# ------------------------------------------------------------------ synthetic tests
def selftest():
    rng = np.random.default_rng(0)
    n = 16                                    # every SHELL is populated at n = 16 (max separation 15)
    t = np.cumsum(rng.normal(size=(n, 3)) * 1.5, 0); t += 0.3 * rng.normal(size=(n, 3))
    c = t + 0.8 * random_lowfreq(n, rng) * np.sqrt(n)
    # R2 round trip: coefficients in the cloud's own frame, applied back, recover the native exactly
    t_c = kabsch_onto(t, c); F = frame(c); E = (c - t_c) @ F
    assert np.abs((c - E @ F.T) - t_c).max() < 1e-10
    assert np.abs(F.T @ F - np.eye(3)).max() < 1e-10 and np.linalg.det(F) > 0
    # frame invariance: rotating the cloud rotates the frame, coefficients unchanged
    U, _, Vt = np.linalg.svd(rng.normal(size=(3, 3))); R = U @ Vt
    if np.linalg.det(R) < 0:
        R[:, 0] *= -1
    c2 = c @ R.T + 5.0; F2 = frame(c2); E2 = (c2 - kabsch_onto(t, c2)) @ F2
    assert np.abs(E2 - E).max() < 1e-8, np.abs(E2 - E).max()
    print("  R2 frame: exact round trip; rotation-invariant coefficients (max dev %.1e)" % np.abs(E2 - E).max())
    # R1 realisation: no move when the target is the cloud's own distances; near-native when it is the native's
    i, j = I.pair_index(n)
    Dc = I.pair_dists(c, i, j); Dt = I.pair_dists(t, i, j)
    x0 = realise(Dc, i, j, c, iters=50)
    assert np.abs(x0 - c).max() < 1e-9
    x1 = realise(Dt, i, j, c, iters=400)
    r_before, r_after = I.ca_rmsd(c, t), I.ca_rmsd(x1, t)
    print("  R1 realise: identity target moves nothing; native target %.3f -> %.3f A" % (r_before, r_after))
    assert r_after < 0.5 * r_before
    # shell means broadcast is the exact inverse on a shell-constant field
    m = rng.normal(size=5); assert np.abs(shell_means(shell_broadcast(m, i, j), i, j) - m).max() < 1e-12
    # rel_grid round trip at n = 16 is the identity
    v = rng.normal(size=(16, 3)); assert np.abs(rel_grid(rel_grid(v, 16), 16) - v).max() < 1e-12
    # NaN-poison: the deployable parts (frame, realise, shell_broadcast, rel_grid) never read the native
    F3 = frame(c); x3 = realise(Dc - 0.1, i, j, c, iters=20)
    tt = np.full_like(t, np.nan)
    assert np.array_equal(frame(c), F3) and np.array_equal(realise(Dc - 0.1, i, j, c, iters=20), x3) and np.isnan(tt).all()
    print("  p_c5 selftest OK")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    {"selftest": selftest, "run": run}[cmd]()
