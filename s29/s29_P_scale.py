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
import math
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
#: THE PRIMARY SET (coordinator, 2026-09-20 00:45, under the CPU-bound one-slot decision):
#: exactly the cells the pre-registered falsifiers F-P1/F-P2/F-P3 need -- PROD, the three
#: candidates, the 8 matched-random draws, CTRL-INV (native-free, so it is in the
#: max-over-K set), CTRL-GLOBAL (likewise) and FLOOR. CTRL-LAM and BOND-LAMFIX are the
#: 2x2 effective-lambda controls, which the prereg reads ONLY if BOND is non-null, so they
#: are DEFERRED with the ORACLE s-grid rather than dropped.
DEFERRED_ARMS = ("CTRL-LAM", "BOND-LAMFIX")


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
    return pr


def shipped_objective(pr, C, cand, pen=None):
    """THE SHIPPED PROJECTION OBJECTIVE, evaluated for an emitted chain against PRODUCTION's own
    cloud at the shipped lam: RMSD(chain, C) + 0.3 * ramah(phi, psi).

    This is what `core.project.fit_multi` minimises over its four generic starts. Evaluating it
    for every arm's emitted chain makes the whole arm set a WIDER MULTI-START of the production
    projection: picking the arm with the lowest value is NATIVE-FREE and strictly improves the
    shipped objective, so it is a deployable operator and not a diagnostic. NATIVE-FREE.
    """
    from core import project as pj
    if pen is None:
        pen = pj.make_penalty("ramah", cand.seq, int(cand.fold))
    phi = np.asarray(pr["phi"], float)[None]
    psi = np.asarray(pr["psi"], float)[None]
    ca = np.asarray(pr["ca"], float)
    return float(I.ca_rmsd(ca, C) + LAM * float(np.asarray(pen(phi, psi)).ravel()[0]))


def run_target(pdb, fac, ref, done, rows_path, floor=True, phase="all"):
    """Every not-yet-done (arm, target) cell for one target, appended as it completes.

    `phase` orders the work, and changes nothing else: "primary" runs the 8 named arms, the 8
    matched-random draws and FLOOR (the cells the pre-registered falsifier needs); "oracle" runs
    the ORACLE s-grid (the ceiling, which is a diagnostic and can land later). Each cell is a
    deterministic function of the cached cloud, its scalar and lam, so the order is immaterial.
    """
    C, cand = production_cloud(pdb)
    assert float(I.ca_rmsd(C, cand.nat_ca)) == float(ref[pdb]["rmsd_cloud"]), \
        f"{pdb}: cached cloud does not reproduce S27's rmsd_cloud exactly"
    A = arms_for(pdb, fac)
    if phase == "primary":
        A = [x for x in A if not x[0].startswith("ORACLE-GRID") and x[0] not in DEFERRED_ARMS]
    elif phase == "oracle":
        A = [x for x in A if x[0].startswith("ORACLE-GRID") or x[0] in DEFERRED_ARMS]
        floor = False
    todo = [(a, s, lm) for (a, s, lm) in A if (pdb, a) not in done]
    if floor and (pdb, "FLOOR") not in done:
        todo.append(("FLOOR", None, LAM))
    n_new = 0
    from core import project as pj
    pen = pj.make_penalty("ramah", cand.seq, int(cand.fold)) if todo else None
    for arm, s, lm in todo:
        t1 = time.time()
        if arm == "FLOOR":
            Cin = C * (1.0 + FLOOR_EPS)
            pr = I.project(Cin, cand.seq, cand.fold, lam=lm)
            s_used = float("nan")
        else:
            Cin = rescale(C, s)
            pr = project_arm(C, cand, s, lm)
            s_used = float(s)
        ca = np.asarray(pr["ca"], float)
        cl = float(I.ca_rmsd(Cin, cand.nat_ca))
        r = dict(pdb=pdb, n=int(cand.n), fold=int(cand.fold), arm=arm, s=s_used, lam=float(lm),
                 #: NATIVE-FREE mechanism quantities: how well the ideal-geometry chain fits the
                 #: cloud it was fitted to, and how far it lands from PRODUCTION's own cloud.
                 fit_resid=float(I.ca_rmsd(ca, Cin)),
                 fit_resid0=float(I.ca_rmsd(ca, C)),
                 #: the SHIPPED objective of this arm's chain against PRODUCTION's cloud:
                 #: native-free, and the quantity a wider multi-start would minimise.
                 obj0=shipped_objective(pr, C, cand, pen),
                 rmsd_chain=float(I.ca_rmsd(ca, cand.nat_ca)),   # ORACLE: post-hoc scoring
                 rmsd_cloud_in=cl,
                 bond_in=mean_bond(Cin), rg_in=rg(Cin),
                 bond_out=mean_bond(ca), rg_out=rg(ca),
                 secs=float(time.time() - t1))
        with open(rows_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(r) + "\n")
        done.add((pdb, arm))
        n_new += 1
    return n_new


def load_rows(path):
    """Read a checkpoint file, tolerating ONE truncated trailing line.

    A shard killed mid-append (re-sharding, a governor kill, a crash) can leave a partial JSON
    object as the last line. Skipping it is correct -- that cell simply was not finished and the
    resume will redo it -- but a partial line anywhere EARLIER means real corruption and must
    raise rather than silently drop completed work.
    """
    out = []
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh]
    lines = [ln for ln in lines if ln]
    for k, line in enumerate(lines):
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            if k == len(lines) - 1:
                print("  note: dropping a truncated trailing line in %s (an unfinished cell; "
                      "it will be recomputed)" % os.path.basename(path), flush=True)
                break
            raise
    return out


def done_keys(paths):
    d = set()
    for p in paths:
        for r in load_rows(p):
            d.add((r["pdb"], r["arm"]))
    return d


def shard_rows_path(shard):
    return os.path.join(RESULTS, f"s29_P_rows_shard{shard}.jsonl")


def cmd_run(shard, nshards, limit=None, pdbs=None, rows_path=None, floor=True, phases=("primary", "oracle")):
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
    #: resume across EVERY shard file, not just this job's own: a consolidated re-run must never
    #: recompute a cell another shard already checkpointed. Each job appends only to its own file.
    allrows = sorted(os.path.join(RESULTS, f) for f in os.listdir(RESULTS)
                     if f.startswith("s29_P_rows_shard") and f.endswith(".jsonl"))
    if rows_path not in allrows:
        allrows.append(rows_path)
    done = done_keys(allrows)
    print("resume: %d (arm, target) cells already checkpointed across %d rows files"
          % (len(done), len(allrows)), flush=True)
    t0 = time.time()
    for ph in phases:
        print(f"== phase {ph} ==", flush=True)
        for k, pdb in enumerate(pdbs):
            n_new = run_target(pdb, fac, ref, done, rows_path, floor=floor, phase=ph)
            print(f"  [{ph} {k+1}/{len(pdbs)}] {pdb} +{n_new} cells "
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


# ============================================================ analysis
def _rand_subset_null(d, k, n_boot=2000, seed_parts=("s29P", "rand18")):
    """The random-k null a stratum claim needs: the percentile of the observed k-subset mean
    among 2000 random k-subsets of the same paired differences."""
    rng = ST._rng(*seed_parts)
    d = np.asarray(d, float)
    n = len(d)
    draws = np.array([d[rng.choice(n, k, replace=False)].mean() for _ in range(n_boot)])
    return draws


def _max_over_k_null(D, n_boot=4000, seed_parts=("s29P", "maxK")):
    """The max-over-K null for the lane's K deployable contrasts, preserving the cross-arm
    correlation: sign-flip the per-target paired differences JOINTLY across arms (a valid
    randomisation under the global null of no effect) and take the max |effect|/MDE."""
    D = np.asarray(D, float)            # (n_targets, K) paired differences arm - PROD
    rng = ST._rng(*seed_parts)
    n, K = D.shape
    out = np.empty(n_boot)
    for t in range(n_boot):
        s = rng.choice([-1.0, 1.0], n)[:, None]
        X = D * s
        eff = X.mean(0)
        mde = ST.MDE_K * X.std(0, ddof=1) / math.sqrt(n)
        out[t] = np.max(np.abs(eff) / mde)
    return out


def cmd_analyse(rows_paths=None, out=None):
    import math as _m
    fac = json.load(open(FACTORS, encoding="utf-8"))
    ref = s27_dis_rows()
    order = fac["order"]
    if rows_paths is None:
        rows_paths = sorted(
            [os.path.join(RESULTS, f) for f in os.listdir(RESULTS)
             if f.startswith("s29_P_rows_shard") and f.endswith(".jsonl")])
    rows = []
    for p in rows_paths:
        rows.extend(load_rows(p))
    cell = {}
    for r in rows:
        cell[(r["pdb"], r["arm"])] = r        # a later row of the same key supersedes (resume)
    arm_names = sorted({a for (_, a) in cell})
    rand_all = [f"CTRL-RAND{d}" for d in range(N_RAND)]
    grid_all = ["PROD" if s == 1.0 else f"ORACLE-GRID{s:.2f}" for s in ORACLE_GRID]
    primary = [a for a in ["PROD", "BOND", "SPAN", "ISO", "CTRL-INV", "CTRL-GLOBAL",
                           "FLOOR"] + rand_all if a in arm_names]
    have = [p for p in order if all((p, a) in cell for a in primary)]
    have_grid = [p for p in order if all((p, a) in cell for a in grid_all)]
    L = []

    def say(s=""):
        L.append(s)
        print(s, flush=True)

    say("LANE P -- THE PROJECTION PRICE. BASIS: BUILT CHAIN (`s12.instrument.project`).")
    say("  rows %d over %d files; %d arms; targets complete on the PRIMARY arms %d / %d; "
        "on the ORACLE grid %d / %d"
        % (len(rows), len(rows_paths), len(arm_names), len(have), len(order),
           len(have_grid), len(order)))
    if len(have) < len(order):
        say("  PRIMARY INCOMPLETE -- missing %s ..."
            % [p for p in order if p not in have][:8])
    pdbs = have
    folds = np.array([cell[(p, "PROD")]["fold"] for p in pdbs])
    ns = np.array([cell[(p, "PROD")]["n"] for p in pdbs])

    def col(arm, key="rmsd_chain"):
        return np.array([cell[(p, arm)].get(key, np.nan) for p in pdbs])

    # ---- the gate, re-asserted on the full run
    prod = col("PROD")
    dref = np.abs(prod - np.array([float(ref[p]["rmsd_chain"]) for p in pdbs]))
    say("  GATE (PROD vs `chain_rows.jsonl :: DIS`): max |diff| %.3e, identical on %d/%d"
        % (dref.max(), int((dref == 0).sum()), len(pdbs)))
    say("")

    rand_cols = rand_all
    grid_cols = grid_all
    gridp = [p for p in have_grid]
    if gridp:
        M_grid = np.column_stack([np.array([cell[(p, a)]["rmsd_chain"] for p in gridp])
                                  for a in grid_cols])
        oracle_scale_g = M_grid.min(1)                             # ORACLE, on `gridp`
        oracle_argmin = np.array([ORACLE_GRID[k] for k in M_grid.argmin(1)])
    else:
        M_grid, oracle_scale_g, oracle_argmin = None, None, None
    #: the ORACLE column aligned to `pdbs` (NaN where the grid is not yet complete)
    gmap = dict(zip(gridp, oracle_scale_g)) if gridp else {}
    oracle_scale = np.array([gmap.get(p, np.nan) for p in pdbs])
    M_rand = np.column_stack([col(a) for a in rand_cols])
    rand_mean = M_rand.mean(1)
    rand_best = M_rand.min(1)                                      # order statistic

    # ---- means table
    say("MEAN BUILT-CHAIN RMSD AND THE PROJECTION PRICE PER ARM")
    say("  price = chain minus THIS ARM's own input cloud; price0 = chain minus PRODUCTION's cloud")
    say("  %-16s %8s %8s %8s %8s %8s %7s %7s %8s %7s"
        % ("arm", "chain", "cloud_in", "price", "price0", "s_mean", "bond_in", "rg_in",
           "bond_out", "rg_out"))
    cloud0 = col("PROD", "rmsd_cloud_in")
    named = [("PROD", prod), ("BOND", col("BOND")), ("SPAN", col("SPAN")), ("ISO", col("ISO")),
             ("CTRL-GLOBAL", col("CTRL-GLOBAL")), ("CTRL-INV", col("CTRL-INV")),
             ("CTRL-LAM", col("CTRL-LAM")), ("BOND-LAMFIX", col("BOND-LAMFIX")),
             ("CTRL-RAND(mean8)", rand_mean), ("CTRL-RAND(best8 ORDER STAT)", rand_best)]
    if np.isfinite(oracle_scale).all():
        named.append(("ORACLE-SCALE", oracle_scale))
    for nm, v in named:
        src = nm if (nm, ) and nm in arm_names else None
        if src:
            say("  %-16s %8.4f %8.4f %+8.4f %+8.4f %8.4f %7.3f %7.3f %8.3f %7.3f"
                % (nm, v.mean(), col(src, "rmsd_cloud_in").mean(),
                   v.mean() - col(src, "rmsd_cloud_in").mean(), v.mean() - cloud0.mean(),
                   col(src, "s").mean(),
                   col(src, "bond_in").mean(), col(src, "rg_in").mean(),
                   col(src, "bond_out").mean(), col(src, "rg_out").mean()))
        else:
            say("  %-16s %8.4f %8s %8s %+8.4f  (derived over columns)"
                % (nm, v.mean(), "-", "-", v.mean() - cloud0.mean()))
    if "FLOOR" in arm_names:
        fl = np.abs(col("FLOOR") - prod)
        say("")
        say("BRANCH-FLIP FLOOR ON THIS CODE PATH (arm FLOOR: C * (1 + 1e-13), same projection)")
        say("  mean |diff| %.4f   max %.4f (%s)   > 0.02 A on %d/%d   > 0.1 A on %d   "
            "mean signed %+.4f"
            % (fl.mean(), fl.max(), pdbs[int(np.argmax(fl))], int((fl > 0.02).sum()), len(fl),
               int((fl > 0.1).sum()), (col("FLOOR") - prod).mean()))
    say("")

    # ---- the contrasts
    C = {}
    say("EVERY ARM AGAINST PRODUCTION (BUILT CHAIN, one code path, one input, n=%d)" % len(pdbs))
    for nm, v in named:
        if nm == "PROD":
            continue
        o = ST.compare(v, prod, folds=folds, names=pdbs,
                       label="P %s - PROD (BUILT CHAIN)%s" % (nm, "  [ORACLE]" if "ORACLE" in nm else ""))
        C[nm] = o
        say(ST.fmt(o)); say("")
    for nm in ("BOND", "SPAN", "ISO", "CTRL-GLOBAL"):
        o = ST.compare(col(nm), rand_mean, folds=folds, names=pdbs,
                       label="P %s - CTRL-RAND mean-of-8 (matched-magnitude derangement of g)" % nm)
        C[nm + "|RAND"] = o
        say(ST.fmt(o)); say("")
    o = ST.compare(col("BOND-LAMFIX"), col("BOND"), folds=folds, names=pdbs,
                   label="P BOND-LAMFIX - BOND (the effective-lambda leg of the 2x2)")
    C["LAMFIX|BOND"] = o
    say(ST.fmt(o)); say("")
    #: THE BRANCH-FLIP REFERENCE. FLOOR is a pure lottery draw: the same cloud perturbed at
    #: 1e-13 relative, i.e. ZERO geometric change, re-projected. If FLOOR - PROD is not itself
    #: zero, the multi-start's branch choice is biased and every arm's effect must be read
    #: against FLOOR rather than against nothing.
    if "FLOOR" in arm_names:
        for nm in ("BOND", "SPAN", "ISO", "CTRL-GLOBAL", "CTRL-INV"):
            if nm not in arm_names:
                continue
            o = ST.compare(col(nm), col("FLOOR"), folds=folds, names=pdbs,
                           label="P %s - FLOOR (against a zero-geometry branch redraw, not "
                                 "against production's particular branch)" % nm)
            C[nm + "|FLOOR"] = o
            say(ST.fmt(o)); say("")

    # ---- THE WIDER MULTI-START (native-free, deployable): every arm's emitted chain is a
    # feasible point of the SHIPPED projection problem (fit production's own cloud C at lam 0.3),
    # so picking the one with the lowest shipped objective obj0 strictly improves production's
    # own optimisation, with no native anywhere. Ties are averaged over the argmin set, never
    # broken by array order.
    ms_arms = [a for a in arm_names if not a.startswith("ORACLE-GRID")]
    OBJ = np.column_stack([col(a, "obj0") for a in ms_arms])
    CH = np.column_stack([col(a) for a in ms_arms])
    #: MS-OBJ / MS-MEAN / MS-ORACLE are DEFERRED with the grid phase (coordinator, 00:45): the
    #: wider multi-start is only reported once its full arm set exists on every analysed target.
    if gridp and len(gridp) == len(pdbs) and np.isfinite(OBJ).all():
        ms_pick, ms_ties = [], []
        for r in range(len(pdbs)):
            v, k = ST.argmin_tied(OBJ[r], CH[r])
            ms_pick.append(v); ms_ties.append(k)
        ms_pick = np.array(ms_pick)
        ms_mean_arm = CH.mean(1)                      # the zero-information selection control
        ms_oracle = CH.min(1)                         # ORACLE selection over the same set
        say("THE WIDER MULTI-START (native-free selection by the SHIPPED projection objective)")
        say("  arms in the set %d; mean tie-set size %.2f; production's own objective is beaten "
            "on %d / %d targets"
            % (len(ms_arms), float(np.mean(ms_ties)),
               int((OBJ.min(1) < col("PROD", "obj0") - 1e-12).sum()), len(pdbs)))
        say("  mean shipped objective: PROD %.4f -> best-over-arms %.4f (a strict improvement "
            "of the quantity production minimises)"
            % (col("PROD", "obj0").mean(), OBJ.min(1).mean()))
        pen_share = col("PROD", "obj0") - col("PROD", "fit_resid0")
        say("  the lam * ramah term contributes %.3e of PROD's objective (mean), max %.3e over "
            "targets: the shipped hinge penalty is inactive at the optimum, so obj0 is the fit "
            "residual to within that -- stated, not assumed"
            % (pen_share.mean(), np.nanmax(np.abs(pen_share))))
        for lab, v in (("MS-OBJ (native-free pick)", ms_pick),
                       ("MS-MEAN (zero-information pick)", ms_mean_arm),
                       ("MS-ORACLE [ORACLE] (pick by RMSD)", ms_oracle)):
            o = ST.compare(v, prod, folds=folds, names=pdbs,
                           label="P %s - PROD (BUILT CHAIN)" % lab)
            C["MS|" + lab.split()[0]] = o
            say(ST.fmt(o)); say("")

    # ---- multiplicity
    #: CTRL-INV (1/g) is native-free, hence DEPLOYABLE, so it enters the multiplicity set
    #: even though the prereg registered it as a direction control. Enlarging K is the
    #: conservative amendment (a higher bar), never the permissive one.
    dep = ("BOND", "SPAN", "ISO", "CTRL-GLOBAL", "CTRL-INV")
    Dcols = [col(a) - prod for a in dep]
    dep_names = list(dep)
    if "MS|MS-OBJ" in C:                      # prereg addendum 4: MS-OBJ is deployable, K = 6
        Dcols.append(ms_pick - prod); dep_names.append("MS-OBJ")
    dep = tuple(dep_names)
    D = np.column_stack(Dcols)
    obs = np.max(np.abs(D.mean(0)) / (ST.MDE_K * D.std(0, ddof=1) / _m.sqrt(len(pdbs))))
    null = _max_over_k_null(D)
    say("MULTIPLICITY: max |effect|/MDE over the %d deployable contrasts %s" % (len(dep), list(dep)))
    say("  observed %.2f   sign-flip max-over-K null p50 %.2f p90 %.2f p95 %.2f   p_value %.3f"
        % (obs, np.percentile(null, 50), np.percentile(null, 90), np.percentile(null, 95),
           float((null >= obs).mean())))
    say("")

    # ---- strata
    F18 = set(I.FAIL18)
    mask = np.array([p in F18 for p in pdbs])
    say("STRATA: FAIL18 (%d here) vs the other %d, with the RANDOM-18 NULL for every claim"
        % (int(mask.sum()), int((~mask).sum())))
    strata_arms = [a for a in ("BOND", "SPAN", "ISO", "ORACLE-SCALE") if a in dict(named)]
    for nm in strata_arms:
        v = dict(named)[nm] if nm == "ORACLE-SCALE" else col(nm)
        d = v - prod
        for lab, mm in (("FAIL18", mask), ("other108", ~mask)):
            o = ST.compare(v[mm], prod[mm], folds=folds[mm], names=[p for p, q in zip(pdbs, mm) if q],
                           label="P %s - PROD, %s%s" % (nm, lab, "  [ORACLE]" if "ORACLE" in nm else ""))
            C["%s|%s" % (nm, lab)] = o
            say("  %s %s effect %+.4f  SE %.4f  MDE %.4f  %+.2fx  fold CI [%+.4f, %+.4f]  "
                "folds %d/%d  %dW/%dL  %s"
                % (nm, lab, o["effect"], o["se"], o["mde"], o["effect_over_mde"],
                   o["ci95_fold"][0], o["ci95_fold"][1], o["folds_same_sign"], o["n_folds"],
                   o["n_better"], o["n_worse"], o["verdict"]))
        draws = _rand_subset_null(d, int(mask.sum()))
        pc = float((draws <= d[mask].mean()).mean())
        say("    random-18 null for %s: observed FAIL18 effect %+.4f sits at percentile %.3f "
            "of 2000 random 18-subsets (null p10/p50/p90 %+.4f/%+.4f/%+.4f)"
            % (nm, d[mask].mean(), pc, np.percentile(draws, 10), np.percentile(draws, 50),
               np.percentile(draws, 90)))
    say("")

    # ---- the ORACLE s-curve and its order-statistic price
    say("THE s-CURVE [ORACLE grid; every number here reads the native and is a CEILING, not a result]")
    bok, cnt = None, {}
    if M_grid is None:
        say("  the grid is not yet complete on any target")
    else:
        say("  n = %d targets complete on the grid" % len(gridp))
        say("  %-8s %8s %8s" % ("s", "mean", "median"))
        for k, (s, a) in enumerate(zip(ORACLE_GRID, grid_cols)):
            v = M_grid[:, k]
            say("  %-8.2f %8.4f %8.4f" % (s, v.mean(), np.median(v)))
        bok = ST.best_of_k_within(M_grid)
        say("  ORACLE-SCALE (per-target argmin over the %d-point grid) mean %.4f; priced: "
            "observed %+.4f, valid null %+.4f (%.0f%%), k_eff %.2f, split-half %+.4f -> %s"
            % (len(ORACLE_GRID), oracle_scale_g.mean(), bok["observed_gain"],
               bok["null_across_targets"], 100 * bok["share_accounted"], bok["k_eff"],
               bok["split_half"], bok["verdict"]))
        cnt = {s: int((oracle_argmin == s).sum()) for s in ORACLE_GRID}
        nties = int(sum((np.isclose(M_grid[r], M_grid[r].min(), rtol=1e-9, atol=1e-12)).sum() > 1
                        for r in range(len(gridp))))
        say("  argmin histogram %s   (the VALUE ORACLE-SCALE reports is the row minimum and is "
            "tie-free; the histogram's ties, %d rows, are broken by grid order and are labelled "
            "so rather than read)" % (cnt, nties))
    bokr = ST.best_of_k_within(M_rand)
    say("  CTRL-RAND best-of-8 priced: observed %+.4f, valid null %+.4f (%.0f%%), split-half "
        "%+.4f -> %s" % (bokr["observed_gain"], bokr["null_across_targets"],
                         100 * bokr["share_accounted"], bokr["split_half"], bokr["verdict"]))
    say("")

    # ---- the contraction check (charter section 16): is a winner just a more compact blob?
    nat_rg, nat_bond = [], []
    for p in pdbs:
        z = np.load(os.path.join(ROOT, "s8", "generate_univ", f"{p}.npz"), allow_pickle=True)
        T = np.asarray(z["nat_ca"], float)          # ORACLE, diagnostic only
        nat_rg.append(rg(T)); nat_bond.append(mean_bond(T))
    nat_rg = np.array(nat_rg); nat_bond = np.array(nat_bond)
    say("CONTRACTION CHECK [ORACLE reference]: native mean Rg %.4f, native mean bond %.4f"
        % (nat_rg.mean(), nat_bond.mean()))
    say("  %-16s %9s %9s %9s" % ("arm", "rg_out", "rg_out/nat", "|log ratio|"))
    for nm, _ in named:
        if nm in arm_names:
            ro = col(nm, "rg_out")
            say("  %-16s %9.4f %9.4f %9.4f"
                % (nm, np.nanmean(ro), np.nanmean(ro / nat_rg),
                   np.nanmean(np.abs(np.log(ro / nat_rg)))))
    say("")

    # ---- the mechanism: does a geometrically consistent cloud admit a better ideal-chain fit?
    say("MECHANISM (NATIVE-FREE): THE PROJECTION RESIDUAL")
    say("  fit_resid  = CA-RMSD(emitted chain, the cloud THAT ARM was fitted to)")
    say("  fit_resid0 = CA-RMSD(emitted chain, PRODUCTION's cloud) -- comparable across arms")
    say("  %-16s %10s %10s" % ("arm", "fit_resid", "fit_resid0"))
    for nm, _ in named:
        if nm in arm_names:
            say("  %-16s %10.4f %10.4f"
                % (nm, np.nanmean(col(nm, "fit_resid")), np.nanmean(col(nm, "fit_resid0"))))
    say("")

    # ---- g diagnostics
    from scipy.stats import spearmanr
    g = np.array([fac["per_target"][p]["g"] for p in pdbs])
    sp = np.array([fac["per_target"][p]["s_span"] for p in pdbs])
    iso = np.array([fac["per_target"][p]["s_iso"] for p in pdbs])
    cloud = col("PROD", "rmsd_cloud_in")
    say("REALISED FACTORS (native-free) AND WHAT THEY TRACK")
    say("  g      mean %.4f sd %.4f  min %.4f max %.4f" % (g.mean(), g.std(ddof=1), g.min(), g.max()))
    say("  s_span mean %.4f sd %.4f ; s_iso mean %.4f sd %.4f"
        % (sp.mean(), sp.std(ddof=1), iso.mean(), iso.std(ddof=1)))
    for nm, x in (("n", ns), ("prod cloud RMSD", cloud), ("prod chain RMSD", prod),
                  ("ORACLE s* [ORACLE]", oracle_argmin)):
        say("  spearman  g vs %-20s %+.3f   s_span %+.3f   s_iso %+.3f"
            % (nm, spearmanr(g, x).statistic, spearmanr(sp, x).statistic, spearmanr(iso, x).statistic))
    say("")

    summ = dict(kind="lane P: the projection price, 126 targets, built chain",
                n=len(pdbs), pdbs=pdbs, arms=arm_names,
                gate_max_diff=float(dref.max()), gate_identical=int((dref == 0).sum()),
                means={nm: float(v.mean()) for nm, v in named},
                price={nm: float(v.mean() - col(nm, "rmsd_cloud_in").mean())
                       for nm, v in named if nm in arm_names},
                n_grid=len(gridp),
                oracle_argmin_hist={str(k): v for k, v in cnt.items()},
                oracle_grid_bok=bok, ctrl_rand_bok=bokr,
                maxk_observed=float(obs), maxk_null_p95=float(np.percentile(null, 95)),
                maxk_p=float((null >= obs).mean()),
                contrasts={k: v for k, v in C.items()},
                floor=(dict(mean=float(np.abs(col("FLOOR") - prod).mean()),
                            max=float(np.abs(col("FLOOR") - prod).max()),
                            n_above_0p02=int((np.abs(col("FLOOR") - prod) > 0.02).sum()))
                       if "FLOOR" in arm_names else None),
                text="\n".join(L))
    ST.save_atomic(out or os.path.join(RESULTS, "s29_P_summary.json"), summ, module_file=__file__)
    print("wrote", out or os.path.join(RESULTS, "s29_P_summary.json"))
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
    r.add_argument("--phase", default="all", choices=["all", "primary", "oracle"])
    an = sub.add_parser("analyse")
    an.add_argument("--out", default=None)
    sub.add_parser("selftest")
    a = ap.parse_args(argv)
    if a.cmd == "factors":
        build_factors(limit=a.limit)
    elif a.cmd == "probe":
        s = cmd_probe()
        sys.exit(0 if s["gate_pass"] else 3)
    elif a.cmd == "run":
        ph = ("primary", "oracle") if a.phase == "all" else (a.phase,)
        cmd_run(a.shard, a.nshards, limit=a.limit, phases=ph)
    elif a.cmd == "analyse":
        cmd_analyse(out=a.out)
    elif a.cmd == "selftest":
        selftest()


if __name__ == "__main__":
    main()
