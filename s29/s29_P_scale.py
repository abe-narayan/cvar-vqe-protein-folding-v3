#!/usr/bin/env python
"""s29/s29_P_scale.py -- LANE P: THE PROJECTION PRICE.

Pre-registration: `s29/PREREG_S29_P.md` (written before any number; arms, falsifier, prior,
the S23 L6 distinction, the probe gate, the multiplicity budget).

THE QUESTION. Production's point cloud is 3.0483 A and its built chain is 3.2126 A: the
projection COSTS +0.159 A. The cloud is geometrically INCONSISTENT with the space it is
projected into -- its mean adjacent CA-CA distance is 2.96 A against the ideal 3.80 -- so the
least-squares fit of a rigid-length chain to it is a BIASED fit. Does removing the geometric
inconsistency BEFORE the projection change the built chain?

NOT S23 L6. L6 closed the RMSD-optimal scale s*, which is read off the NATIVE and is a property
of the (pool, reference-structure) pair. Every deployable scale here is read off the CLOUD's own
bonds or off the leave-fold-out distogram posterior, and the emitted object is an ideal-geometry
chain whose CA-CA distance is 3.80 by construction -- so a rescale of the INPUT cannot rescale
the OUTPUT, only change which ideal chain is nearest. See PREREG section 2.

ORACLE: `oracle_scale` and every `rmsd_chain` read the native. The native is read ONLY to score
an emitted structure post hoc, and by the ORACLE-SCALE arm, which is labelled ORACLE in every
line it appears. No native quantity chooses a deployable parameter.

    python s29/s29_P_scale.py factors  [--limit N]      # cloud cache + every scale factor
    python s29/s29_P_scale.py probe                     # the 6 registered targets, 1e-9 gate
    python s29/s29_P_scale.py run --shard K --nshards N # the 126, per-(arm, target) checkpoint
    python s29/s29_P_scale.py analyse                   # concatenate the shards, ST.fmt blocks
    python s29/s29_P_scale.py selftest
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s24 import d_harness as H             # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import run_vqe_chain as RV        # noqa: E402

RESULTS = os.path.join(HERE, "results")
CLOUDS = os.path.join(RESULTS, "s29_P_clouds")
FACTORS = os.path.join(RESULTS, "s29_P_factors.json")
S27_CHAIN = os.path.join(ROOT, "s27", "results", "chain_rows.jsonl")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(CLOUDS, exist_ok=True)

#: the ideal virtual CA-CA distance of the builder the projection emits into
#: (`core.project.build_ca_exact`; s23 L1 measured the native panel at 3.812).
IDEAL_BOND = 3.80
LAM = 0.3
N_RAND = 8
#: the ORACLE grid (PREREG section 5). 1.00 is PROD and is not re-projected.
ORACLE_GRID = (0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30)
PROBE = ("1A13", "1CS9", "2LNG", "2NB7", "9BFL", "9BAF")
#: the relative perturbation used by arm FLOOR to measure the branch-flip floor on THIS path
FLOOR_EPS = 1e-13


# ============================================================ pinned reference rows
def s27_dis_rows():
    """S27's production rows, config DIS: the anchor both the cloud and the chain gate against."""
    ref = {}
    with open(S27_CHAIN, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if r.get("config") == "DIS":
                ref[r["pdb"]] = r
    return ref


def targets():
    return sorted(t["pdb"] for t in I.targets())


# ============================================================ the production cloud
def production_cloud(pdb, cache=True):
    """THE DEPLOYED POINT CLOUD, reproduced on `s27/s28_B_prodcheck.py :: project_production`'s
    code path (bit-exact on 126/126, `s27/results/s28_B_prodcheck.json`).

    Returns (C, cand). C is cached as float64 npz -- lossless, so the cached array is the
    array the deployed readout emitted.
    """
    f = os.path.join(CLOUDS, f"{pdb}.npz")
    cand, ch, _ = RP.channels_for(pdb)
    if cache and os.path.exists(f):
        z = np.load(f)
        return np.asarray(z["C"], float), cand
    key = RP.rng_for(pdb, "tiekey").random(cand.k)
    E = RV.energy_for("DIS", ch, pdb, key)
    top = np.lexsort((key, E))[:RP.M]
    C, _ = H.readout_uniform(cand, top)
    C = np.asarray(C, float)
    if cache:
        tmp = f + f".tmp{os.getpid()}.npz"
        np.savez(tmp, C=C)
        os.replace(tmp, f)
    return C, cand


# ============================================================ geometry
def mean_bond(C):
    """Mean adjacent (virtual) CA-CA distance of a CA trace."""
    C = np.asarray(C, float)
    return float(np.linalg.norm(C[1:] - C[:-1], axis=1).mean())


def rg(C):
    """Radius of gyration from coordinates."""
    C = np.asarray(C, float)
    return float(np.sqrt(((C - C.mean(0)) ** 2).sum(1).mean()))


def rg_from_pairs(D, n):
    """The classical identity Rg^2 = (1/2N^2) sum_ij d_ij^2, with D the FULL ordered pair
    distance matrix (n, n).  Equals `rg` exactly for an embeddable D."""
    D = np.asarray(D, float)
    return float(np.sqrt((D ** 2).sum() / (2.0 * n * n)))


def rescale(C, s):
    """C_s = mu + s (C - mu).  s = 1.0 is the exact float identity."""
    C = np.asarray(C, float)
    if s == 1.0:
        return C
    mu = C.mean(0)
    return mu + s * (C - mu)


# ============================================================ the posterior median map
def posterior_median(dg):
    """The L1-Bayes (median) distance per pair: argmin over the grid of the stored L1 risk.

    `dg['risk']` is the per-pair L1 Bayes risk on `dg['grid']` (2..40 by 0.05), weighted by a
    positive per-pair weight, so the argmin is the posterior median regardless of the weight.
    Returns (i, j, m) with |i - j| >= 2 (the pair set the distogram supplies).
    """
    grid = np.asarray(dg["grid"], float)
    risk = np.asarray(dg["risk"], float)
    m = grid[np.argmin(risk, axis=1)]
    return np.asarray(dg["i"], int), np.asarray(dg["j"], int), np.asarray(m, float)


def cloud_pair_dists(C, i, j):
    C = np.asarray(C, float)
    return np.linalg.norm(C[i] - C[j], axis=1)


# ============================================================ the deployable scales
def g_bond(C):
    """BOND: the factor that makes the cloud's mean adjacent CA-CA distance the ideal 3.80 A.

    NATIVE-FREE by construction: it reads the cloud and one covalent constant.
    """
    b = mean_bond(C)
    return float(IDEAL_BOND / b)


def s_span(C, dg):
    """SPAN: the factor that makes the cloud's Rg equal the Rg the posterior's own median map
    implies, through Rg^2 = (1/2N^2) sum_ij d_ij^2 with the |i-j| = 1 terms completed at the
    ideal 3.80 A (the posterior supplies |i-j| >= 2 only).  NATIVE-FREE.
    """
    C = np.asarray(C, float)
    n = len(C)
    i, j, m = posterior_median(dg)
    D = np.zeros((n, n), float)
    D[i, j] = m
    D[j, i] = m
    k = np.arange(n - 1)
    D[k, k + 1] = IDEAL_BOND
    D[k + 1, k] = IDEAL_BOND
    rg_post = rg_from_pairs(D, n)
    return float(rg_post / rg(C)), float(rg_post)


def s_spanr(C, dg):
    """The restricted SPAN diagnostic: the RMS ratio over the posterior's OWN pairs, with no
    3.80 injected.  Reported, never projected (PREREG section 5).  NATIVE-FREE.
    """
    i, j, m = posterior_median(dg)
    d = cloud_pair_dists(C, i, j)
    return float(np.sqrt((m ** 2).sum() / (d ** 2).sum()))


def s_iso(C, dg):
    """ISO: the one-parameter least-squares calibration of the cloud's distance map onto the
    posterior median map, argmin_s sum (s d_ij - m_ij)^2 = <d, m> / <d, d>, over the
    posterior's own pair set.  NATIVE-FREE.
    """
    i, j, m = posterior_median(dg)
    d = cloud_pair_dists(C, i, j)
    return float(float(d @ m) / float(d @ d))


def iso_residual(C, dg):
    """The least-squares residual ISO leaves, at its own optimal scalar: the native-free
    discrepancy that would become a cost function if ISO produced an endpoint positive
    (PREREG section 9)."""
    i, j, m = posterior_median(dg)
    d = cloud_pair_dists(C, i, j)
    s = float(d @ m) / float(d @ d)
    return float(np.sqrt(((s * d - m) ** 2).mean()))


def cost_isoresid(W, ctx):
    """Lane D meter contract: f(W, ctx) -> (m,) float, lower is better.  The ISO residual of
    each cloud against the posterior median map.  Reads only ctx.dg and the clouds; never
    ctx.nat_ca (NaN by construction in the meter)."""
    W = np.asarray(W, float)
    if W.ndim == 2:
        W = W[None]
    dg = ctx.dg if hasattr(ctx, "dg") else ctx["dg"]
    i, j, m = posterior_median(dg)
    out = np.empty(len(W))
    for a in range(len(W)):
        d = np.linalg.norm(W[a][i] - W[a][j], axis=1)
        den = float(d @ d)
        s = (float(d @ m) / den) if den > 0 else 1.0
        out[a] = np.sqrt(((s * d - m) ** 2).mean())
    return out


# ============================================================ the controls
def derangements(gs, n_draw=N_RAND):
    """CTRL-RAND: `n_draw` derangements of the realised g across the 126 targets.  The multiset
    of factors is EXACTLY the multiset of realised g (matched magnitude), and no target keeps
    its own (zero per-target information).  Seeded from sha256, never from array order.
    """
    gs = np.asarray(gs, float)
    n = len(gs)
    out = []
    for dr in range(n_draw):
        h = hashlib.sha256(f"s29P|ctrlrand|{dr}".encode()).digest()
        rng = np.random.default_rng(int.from_bytes(h[:8], "little"))
        for _ in range(10000):
            p = rng.permutation(n)
            if not np.any(p == np.arange(n)):
                break
        else:
            raise RuntimeError("no derangement found")
        out.append(p)
    return np.asarray(out, int)


# ============================================================ the arm table
def arms_for(pdb, fac):
    """Every (arm, s, lam) this lane projects for one target, from the pre-computed factors."""
    f = fac["per_target"][pdb]
    gmean = fac["g_mean"]
    A = [("PROD", 1.0, LAM),
         ("BOND", f["g"], LAM),
         ("SPAN", f["s_span"], LAM),
         ("ISO", f["s_iso"], LAM),
         ("CTRL-INV", 1.0 / f["g"], LAM),
         ("CTRL-GLOBAL", gmean, LAM),
         ("CTRL-LAM", 1.0, LAM * gmean),
         ("BOND-LAMFIX", f["g"], LAM * f["g"])]
    for d in range(N_RAND):
        A.append((f"CTRL-RAND{d}", f["rand"][d], LAM))
    for s in ORACLE_GRID:
        if s == 1.0:
            continue
        A.append((f"ORACLE-GRID{s:.2f}", s, LAM))
    return A


DEPLOYABLE = ("PROD", "BOND", "SPAN", "ISO", "CTRL-GLOBAL", "CTRL-LAM", "BOND-LAMFIX")


# ============================================================ factors stage
def build_factors(limit=None, verbose=True):
    """Cache the 126 production clouds, assert each reproduces S27's `rmsd_cloud` EXACTLY, and
    compute every scale factor.  Writes `s29/results/s29_P_factors.json`."""
    ref = s27_dis_rows()
    pdbs = targets()
    if limit:
        pdbs = pdbs[:limit]
    per, bad = {}, []
    t0 = time.time()
    for k, pdb in enumerate(pdbs):
        C, cand = production_cloud(pdb)
        rc = float(I.ca_rmsd(C, cand.nat_ca))          # ORACLE: post-hoc scoring of the cloud
        d = abs(rc - float(ref[pdb]["rmsd_cloud"]))
        if d != 0.0:
            bad.append((pdb, d))
        dg = I.distogram(pdb, cand.seq, cand.fold)
        sp, rgp = s_span(C, dg)
        per[pdb] = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold),
                        g=g_bond(C), s_span=sp, s_iso=s_iso(C, dg), s_spanr=s_spanr(C, dg),
                        rg_post=rgp, rg_cloud=rg(C), bond_cloud=mean_bond(C),
                        iso_resid=iso_residual(C, dg),
                        rmsd_cloud=rc, ref_cloud=float(ref[pdb]["rmsd_cloud"]),
                        d_cloud=d, ref_chain=float(ref[pdb]["rmsd_chain"]))
        if verbose and (k + 1) % 10 == 0:
            print(f"  [{k+1}/{len(pdbs)}] {pdb} g {per[pdb]['g']:.4f} "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)
    order = [p for p in pdbs]
    gs = np.array([per[p]["g"] for p in order])
    P = derangements(gs)
    for a, p in enumerate(order):
        per[p]["rand"] = [float(gs[P[d][a]]) for d in range(N_RAND)]
    out = dict(kind="lane P scale factors; clouds reproduced on s28_B_prodcheck's code path",
               n=len(order), order=order, g_mean=float(gs.mean()), g_sd=float(gs.std(ddof=1)),
               ideal_bond=IDEAL_BOND, lam=LAM, oracle_grid=list(ORACLE_GRID),
               n_cloud_mismatch=len(bad), cloud_mismatch=bad, per_target=per)
    ST.save_atomic(FACTORS, out, module_file=__file__)
    if verbose:
        print(f"g: mean {out['g_mean']:.4f} sd {out['g_sd']:.4f} "
              f"min {gs.min():.4f} max {gs.max():.4f}; cloud mismatches {len(bad)}")
        print("wrote", FACTORS)
    return out


# ============================================================ projection stage
def project_arm(C, cand, s, lam):
    """ONE projection.  The only place this lane calls the pipeline's stage 3b."""
    pr = I.project(rescale(C, s), cand.seq, cand.fold, lam=lam)
    ca = np.asarray(pr["ca"], float)
    return ca


def run_target(pdb, fac, ref, done, rows_path, floor=True):
    """Every not-yet-done (arm, target) cell for one target, appended as it completes."""
    C, cand = production_cloud(pdb)
    assert float(I.ca_rmsd(C, cand.nat_ca)) == float(ref[pdb]["rmsd_cloud"]), \
        f"{pdb}: cached cloud does not reproduce S27's rmsd_cloud exactly"
    todo = [(a, s, lm) for (a, s, lm) in arms_for(pdb, fac) if (pdb, a) not in done]
    if floor and (pdb, "FLOOR") not in done:
        todo.append(("FLOOR", None, LAM))
    n_new = 0
    for arm, s, lm in todo:
        t1 = time.time()
        if arm == "FLOOR":
            Cp = C * (1.0 + FLOOR_EPS)
            ca = I.project(Cp, cand.seq, cand.fold, lam=lm)["ca"]
            cl = float(I.ca_rmsd(Cp, cand.nat_ca))
            s_used = float("nan")
        else:
            Cs = rescale(C, s)
            ca = project_arm(C, cand, s, lm)
            cl = float(I.ca_rmsd(Cs, cand.nat_ca))
            s_used = float(s)
        ca = np.asarray(ca, float)
        r = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), arm=arm, s=s_used, lam=float(lm),
                 rmsd_chain=float(I.ca_rmsd(ca, cand.nat_ca)),   # ORACLE: post-hoc scoring
                 rmsd_cloud_in=cl,
                 bond_in=mean_bond(rescale(C, s) if arm != "FLOOR" else C),
                 rg_in=rg(rescale(C, s) if arm != "FLOOR" else C),
                 bond_out=mean_bond(ca), rg_out=rg(ca),
                 secs=float(time.time() - t1))
        with open(rows_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(r) + "\n")
        done.add((pdb, arm))
        n_new += 1
    return n_new


def load_rows(path):
    out = []
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def done_keys(paths):
    d = set()
    for p in paths:
        for r in load_rows(p):
            d.add((r["pdb"], r["arm"]))
    return d


def shard_rows_path(shard):
    return os.path.join(RESULTS, f"s29_P_rows_shard{shard}.jsonl")


def cmd_run(shard, nshards, limit=None, pdbs=None, rows_path=None, floor=True):
    fac = json.load(open(FACTORS, encoding="utf-8"))
    ref = s27_dis_rows()
    if pdbs is None:
        allp = fac["order"]
        if limit:
            allp = allp[:limit]
        if nshards > 1:
            # contiguous quarters of the pinned pdb-sorted order
            q = int(np.ceil(len(allp) / nshards))
            pdbs = allp[shard * q:(shard + 1) * q]
        else:
            pdbs = allp
    rows_path = rows_path or shard_rows_path(shard)
    done = done_keys([rows_path])
    t0 = time.time()
    for k, pdb in enumerate(pdbs):
        n_new = run_target(pdb, fac, ref, done, rows_path, floor=floor)
        print(f"  [{k+1}/{len(pdbs)}] {pdb} +{n_new} cells "
              f"({(time.time()-t0)/60:.1f} min elapsed)", flush=True)
    print(f"shard {shard}: {len(pdbs)} targets, rows {rows_path}, {(time.time()-t0)/60:.1f} min")


# ============================================================ the probe and its gate
def cmd_probe():
    """The 6 registered targets; GATE: PROD reproduces S27's built chain to < 1e-9."""
    if not os.path.exists(FACTORS):
        build_factors()
    fac = json.load(open(FACTORS, encoding="utf-8"))
    ref = s27_dis_rows()
    rows_path = os.path.join(RESULTS, "s29_P_rows_probe.jsonl")
    done = done_keys([rows_path])
    t0 = time.time()
    for k, pdb in enumerate(PROBE):
        n_new = run_target(pdb, fac, ref, done, rows_path)
        print(f"  [{k+1}/6] {pdb} +{n_new} cells ({(time.time()-t0)/60:.1f} min)", flush=True)
    rows = load_rows(rows_path)
    prod = {r["pdb"]: r["rmsd_chain"] for r in rows if r["arm"] == "PROD"}
    gate = {p: abs(prod[p] - float(ref[p]["rmsd_chain"])) for p in PROBE}
    mx = max(gate.values())
    summ = dict(kind="lane P 6-target probe; GATE = PROD reproduces chain_rows.jsonl::DIS to <1e-9",
                targets=list(PROBE), gate_diff=gate, gate_max=float(mx),
                gate_pass=bool(mx < 1e-9), n_identical=int(sum(v == 0.0 for v in gate.values())),
                rows=rows,
                secs_per_projection=float(np.mean([r["secs"] for r in rows])),
                n_cells=len(rows))
    by_arm = {}
    for r in rows:
        by_arm.setdefault(r["arm"], []).append(r["rmsd_chain"])
    summ["probe_means"] = {a: float(np.mean(v)) for a, v in sorted(by_arm.items())}
    ST.save_atomic(os.path.join(RESULTS, "s29_P_probe.json"), summ, module_file=__file__)
    print(json.dumps({k: v for k, v in summ.items() if k != "rows"}, indent=1))
    print("GATE", "PASS" if summ["gate_pass"] else "FAIL", f"max |diff| {mx:.3e}")
    return summ


# ============================================================ selftest
def selftest():
    rng = np.random.default_rng(0)
    C = rng.normal(size=(12, 3)) * 3.0
    assert rescale(C, 1.0) is C or np.array_equal(rescale(C, 1.0), C)
    n = len(C)
    D = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=-1)
    assert abs(rg_from_pairs(D, n) - rg(C)) < 1e-12, "Rg identity"
    g = g_bond(C)
    assert abs(mean_bond(rescale(C, g)) - IDEAL_BOND) < 1e-12, "BOND"
    gs = np.linspace(1.0, 1.5, 20)
    P = derangements(gs, 3)
    for p in P:
        assert not np.any(p == np.arange(len(gs)))
        assert sorted(gs[p].tolist()) == sorted(gs.tolist())
    print("selftest ok")


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("factors"); f.add_argument("--limit", type=int, default=None)
    sub.add_parser("probe")
    r = sub.add_parser("run")
    r.add_argument("--shard", type=int, default=0)
    r.add_argument("--nshards", type=int, default=1)
    r.add_argument("--limit", type=int, default=None)
    sub.add_parser("selftest")
    a = ap.parse_args(argv)
    if a.cmd == "factors":
        build_factors(limit=a.limit)
    elif a.cmd == "probe":
        s = cmd_probe()
        sys.exit(0 if s["gate_pass"] else 3)
    elif a.cmd == "run":
        cmd_run(a.shard, a.nshards, limit=a.limit)
    elif a.cmd == "selftest":
        selftest()


if __name__ == "__main__":
    main()
