"""s19/agentC_lib.py -- AGENT C (LEGACY / AMBER / PHYSICS), Sprint 19 shared machinery.

Three things live here and nothing else:

  1. THE COMMON FRAME of the deployed averaging operator, and the Krogh-Vedelsby ambiguity
     decomposition on it.  This is the instrument for question 1.
  2. The gates (score gates, matched-random gates, diversity-preserving gates) -- all
     native-free, all operating on the SAME 75 windows.
  3. Artefact / statistics plumbing, imported from `s18.phys_lib` so Sprint 18 and Sprint 19
     numbers are produced by literally the same estimator.

--------------------------------------------------------------------------------------
THE DECOMPOSITION, AND WHY IT IS AN IDENTITY AND NOT A DISCOVERY  (label: EXACT)

The deployed operator emits  Ybar = (1/m) sum_i Y_i  where Y_i is candidate i's ideal-geometry
backbone rebuild, moved into ONE COMMON FRAME (the medoid frame of the surviving window set --
`s15.phys_repl.averaged_backbone_from`).  Let T be the native CA trace optimally superposed
onto Ybar.  Writing e_i = Y_i - T and b = Ybar - T = mean_i e_i, the parallel-axis theorem gives,
exactly, with ||.||^2 the per-residue mean square:

    readout^2  =  ||b||^2  =  (1/m) sum_i ||e_i||^2   -   (1/m) sum_i ||Y_i - Ybar||^2
                              [------ E_mem^2 ------]     [----- D^2  (ambiguity) -----]

    readout^2 = E_mem^2 - D^2.                                                      (KV)

This is a theorem.  What is NOT a theorem, and is the whole question, is HOW A GATE MOVES THE
TWO RIGHT-HAND TERMS.  A gate that raises set quality (lowers E_mem) can still raise the readout
if it destroys D by more.

THE TRAP THIS MODULE EXISTS TO AVOID.  A prior sprint broke this identity by 18% by mixing
FREE-SUPERPOSITION member RMSDs (`I.kabsch_rmsd_batch(W, nat)`, each member optimally rotated
onto the native on its own) with COMMON-FRAME ensemble quantities.  Those are different numbers:

    d_i        = free-superposition CA-RMSD of member i to the native          (what the ledger's
                 operator law d_out = 1.16 d_set_mean + 0.04 d_set_best consumes)
    sqrt(mm_i) = COMMON-FRAME error of member i to the native                  (what KV consumes)

and mm_i >= d_i^2 always, because d_i is the minimum over rigid motions -- PROVIDED both are
measured on the SAME object.  They are not, by default: the operator law's `d` is the RAW WINDOW
CA trace `W[i]`, while the operator averages the IDEAL-GEOMETRY REBUILD of (phi_i, psi_i).  Those
differ by the project's 0.347 A rebuild floor, and using the raw-window `d` in (KV3) makes
FRAME^2 come out NEGATIVE (measured: -0.059 A^2 in gate GC0's first draft).  This module
therefore carries BOTH, never one in the other's place:

    d_reb_i    free superposition of the REBUILD Y_i onto the native  -> enters (KV3)
    d_win_i    free superposition of the WINDOW W[i] onto the native  -> the operator law's input

The gap

    FRAME^2 := (1/m) sum_i mm_i  -  (1/m) sum_i d_reb_i^2   >= 0

is the ensemble's internal misalignment: how much of each member's error is a rigid mismatch that
the common frame cannot remove because the frame is shared.  So the FULL accounting is

    readout^2  =  <d_reb^2>  +  FRAME^2  -  D^2                                      (KV3)

with all three terms measured, and only the first visible to the operator law.  `kv_of` returns
all of them and asserts (KV) to 1e-9 A^2 on every call.

--------------------------------------------------------------------------------------
NATIVE INFORMATION.  T (hence readout, E_mem, FRAME) is ORACLE and appears only in evaluation and
in labelled ORACLE diagnostics.  D and every gate in this module are native-free.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                   # noqa: E402
from s14.avgspace import top75_windows                            # noqa: E402
from s14.ener_avgrefine import _kabsch                            # noqa: E402
from s15 import seed as SD                                        # noqa: E402
from s16 import energy_lib as EL                                  # noqa: E402
from s18 import phys_lib as PL                                    # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

ATOMS = ("N", "CA", "C", "O", "CB")

#: The eleven genuine Legacy components (`core.energy`, DEFAULT_WEIGHTS, never fitted).
LEG_TERMS = EL.LEG_TERMS


# ==========================================================================
# 1.  THE COMMON FRAME OF THE DEPLOYED OPERATOR
# ==========================================================================
def common_frame_members(W, PHI, PSI):
    """The deployed operator's members and mean, in ITS OWN common frame.

    Line-for-line `s15.phys_repl.averaged_backbone_from`, except that the per-member moved
    coordinates are RETURNED instead of being consumed by `.mean(0)`.  Verified against it in
    gate GC0: `Y.mean(0)` reproduces `avg["CA"]` to machine precision, so this is the same
    operator, not a re-implementation of it.

    Returns (Y, Ybar, C_ca, dev) with Y (m, n, 3) the members' CA in the common frame.
    """
    from core import geometry as geo
    W = np.asarray(W, float)
    C_ca, b = I.coordinate_average(W)
    Wm = I.superpose_batch(W, W[b])
    full = geo.build_backbone_batch(np.asarray(PHI, float), np.asarray(PSI, float))
    bca = np.asarray(full["CA"], float)
    R, Am, Bm = _kabsch(bca, Wm)
    Y = np.einsum("bij,bnj->bni", R, bca - Am) + Bm
    Ybar = Y.mean(0)
    return Y, Ybar, C_ca, float(np.abs(Ybar - C_ca).max())


def kv_of(Y, nat, d_win=None):
    """The Krogh-Vedelsby decomposition of one survivor set.  (KV) is an identity: asserted.

    `Y` (m, n, 3) common-frame members; `nat` (n, 3) the native CA trace (ORACLE);
    `d_win` optional (m,) free-superposition CA-RMSDs of the RAW WINDOWS (the operator law's
    input), carried beside -- never inside -- the identity.  `d_reb` (the free-superposition
    RMSD of the SAME rebuilt members that enter the average) is computed here, because only
    that one makes FRAME^2 >= 0.
    """
    Y = np.asarray(Y, float)
    m, n, _ = Y.shape
    Ybar = Y.mean(0)
    #: the native in the readout's own frame -- the metric superposes the emitted mean onto the
    #: native, so the same rigid motion defines every member's common-frame error.
    T = I.superpose_batch(np.asarray(nat, float)[None], Ybar)[0]
    e = Y - T                                       # (m, n, 3)
    readout2 = float((np.sum((Ybar - T) ** 2)) / n)
    emem2 = float(np.sum(e ** 2) / (m * n))
    D2 = float(np.sum((Y - Ybar) ** 2) / (m * n))
    resid = abs(readout2 - (emem2 - D2))
    if resid > 1e-9 * max(1.0, emem2):
        raise AssertionError(f"KV identity broken by {resid:.3e} A^2")
    out = {"m": int(m),
           "readout": float(np.sqrt(max(readout2, 0.0))),
           "E_mem": float(np.sqrt(max(emem2, 0.0))),
           "D": float(np.sqrt(max(D2, 0.0))),
           "readout2": readout2, "E_mem2": emem2, "D2": D2,
           "kv_resid": float(resid)}
    #: (KV3): the free-superposition RMSD of the SAME rebuilt members that enter the average.
    dreb = np.asarray(I.kabsch_rmsd_batch(Y, np.asarray(nat, float)), float)
    out["d_reb_mean"] = float(dreb.mean())
    out["d_reb_best"] = float(dreb.min())
    out["d_reb2_mean"] = float((dreb ** 2).mean())
    out["FRAME2"] = float(emem2 - (dreb ** 2).mean())
    out["FRAME"] = float(np.sqrt(max(out["FRAME2"], 0.0)))
    if d_win is not None:
        df = np.asarray(d_win, float)
        out["d_mean"] = float(df.mean())
        out["d_best"] = float(df.min())
        out["d2_mean"] = float((df ** 2).mean())
    #: mean pairwise cosine of the member ERROR vectors -- the interpretable face of D/E_mem.
    #: (identity: mean_{i != j} <e_i,e_j> = (m ||b||^2 - E_mem^2 * ... ) ; reported as a
    #: diagnostic, the decisive quantities are readout2/E_mem2/D2 above.)
    E = e.reshape(m, -1)
    nrm = np.linalg.norm(E, axis=1)
    ok = nrm > 1e-12
    if ok.sum() >= 2:
        U = E[ok] / nrm[ok][:, None]
        G = U @ U.T
        iu = np.triu_indices(len(U), 1)
        out["err_cos"] = float(G[iu].mean())
    else:
        out["err_cos"] = float("nan")
    out["shared_frac"] = float(readout2 / emem2) if emem2 > 0 else float("nan")
    return out


# ==========================================================================
# 2.  GATES.  Every gate is native-free and returns indices into the SAME 75 windows.
# ==========================================================================
def keep_lowest(score, n_keep):
    """Indices of the `n_keep` LOWEST-scoring candidates.  Identical to `s18.phys_down._keep`,
    including its deterministic index tie-break, so the Sprint-18 arms are reproduced exactly."""
    s = np.asarray(score, float)
    s = np.where(np.isfinite(s), s, np.inf)
    return np.sort(np.lexsort((np.arange(len(s)), s))[:n_keep])


def gate_random(K, n_keep, rng):
    return np.sort(rng.permutation(K)[:n_keep])


def _farthest_point(P, n_keep, start):
    """Greedy max-min selection on a precomputed pairwise-distance matrix.  Deterministic."""
    P = np.asarray(P, float)
    sel = [int(start)]
    dmin = P[start].copy()
    while len(sel) < n_keep:
        dmin[sel] = -np.inf
        j = int(np.argmax(dmin))
        sel.append(j)
        dmin = np.minimum(dmin, P[j])
    return np.sort(np.asarray(sel, int))


def gate_score_then_spread(score, P, n_keep, pre=2.0):
    """DIVERSITY-PRESERVING SCORE GATE.  Native-free.

    Reject on the score down to `pre * n_keep` survivors, then pick `n_keep` of those by greedy
    farthest-point on CA-RMSD, seeded at the best-scoring survivor.  This is a rejection rule of
    exactly the same COUNT as the plain score gate, using exactly the same score, that is built
    not to collapse the ensemble's spread.
    """
    K = len(np.asarray(score))
    n_pre = min(K, max(n_keep, int(round(pre * n_keep))))
    pre_idx = keep_lowest(score, n_pre)
    sub = _farthest_point(np.asarray(P, float)[np.ix_(pre_idx, pre_idx)], n_keep,
                          int(np.argmin(np.asarray(score, float)[pre_idx])))
    return np.sort(pre_idx[sub])


def gate_cluster_best(score, P, n_keep, rng):
    """DIVERSITY-PRESERVING SCORE GATE, second flavour.  Native-free.

    Partition the 75 into `n_keep` clusters by CA-RMSD (farthest-point seeding + one Lloyd-style
    reassignment) and keep each cluster's BEST-scoring member.  The score does all the rejecting
    WITHIN a structural neighbourhood; the ensemble keeps one representative per neighbourhood.
    """
    P = np.asarray(P, float)
    K = len(P)
    seeds = _farthest_point(P, n_keep, int(rng.integers(K)))
    lab = np.argmin(P[np.ix_(seeds, np.arange(K))], axis=0)
    s = np.asarray(score, float)
    out = []
    for c in range(n_keep):
        mem = np.flatnonzero(lab == c)
        if len(mem) == 0:
            continue
        out.append(int(mem[np.lexsort((mem, s[mem]))[0]]))
    out = sorted(set(out))
    #: a collapsed cluster can leave the set short; top it up by score, never at random, so the
    #: arm stays deterministic and its COUNT stays matched.
    if len(out) < n_keep:
        for j in keep_lowest(s, K):
            if int(j) not in out:
                out.append(int(j))
            if len(out) == n_keep:
                break
    return np.sort(np.asarray(out[:n_keep], int))


# ==========================================================================
# 3.  SCORES on the shipped top-75, all native-free
# ==========================================================================
def legacy_scores(seq, PHI, PSI):
    comp = EL.legacy_components_of_windows(seq, PHI, PSI)
    tot = np.asarray(EL.legacy_total_from(comp), float)
    return comp, tot


def helix_score(W, n):
    """ZERO-INFORMATION reference score: CA-RMSD of each candidate to a constant ideal helix."""
    from s17 import phys_lib as P17
    hb = P17.helix_backbone(int(n))
    return np.asarray(I.kabsch_rmsd_batch(np.asarray(W, float),
                                          np.asarray(hb["CA"], float)), float)


# ==========================================================================
# 4.  GATES (verification), run before any scientific number is quoted
# ==========================================================================
def concentration(d, n_drop=5, n_sim=4000, seed=0):
    """Is a mean carried by a few targets?  Drop-top compared to a UNIFORM-EFFECT null.

    The recorded trap (`median-vs-mean-is-the-free-warning`): a raw drop-top threshold is NOT a
    valid test, because dropping the largest `k` of ANY noisy sample shrinks its mean.  The null
    here is a uniform effect of the observed size plus resampled centred noise, so the reported
    percentile is the percentile of the OBSERVED drop-top mean within what a genuinely uniform
    effect would produce.  A low percentile means the effect is concentrated.
    """
    d = np.asarray(d, float)
    d = d[np.isfinite(d)]
    n = len(d)
    mu = float(d.mean())
    obs = float(np.sort(d)[:max(n - n_drop, 1)].mean())      # drop the `n_drop` LARGEST
    rng = np.random.default_rng(seed)
    e = d - mu
    sim = np.empty(n_sim)
    for b in range(n_sim):
        x = mu + e[rng.integers(0, n, n)]
        sim[b] = np.sort(x)[:max(n - n_drop, 1)].mean()
    pct = float((sim <= obs).mean())
    return {"n": int(n), "mean": mu, "median": float(np.median(d)), "sd": float(d.std(ddof=1)),
            "W": int((d < 0).sum()), "L": int((d > 0).sum()),
            "drop_top_mean": obs, "null_pct": pct,
            "verdict": ("CONCENTRATED" if pct < 0.05 else
                        "not concentrated" if pct > 0.20 else "borderline")}


def gate_GC0(n_targets=6, verbose=True):
    """GC0 -- `common_frame_members` IS the deployed operator, and (KV) holds to 1e-12.

    Checks, on the full 75 and on a gated 37:
      a) Y.mean(0) == averaged_backbone_from(...)[0]["CA"]        (max |delta|)
      b) ca_rmsd(Y.mean(0), nat) == sqrt(readout2)                (the metric IS the readout)
      c) readout2 == E_mem2 - D2                                  (the identity)
      d) FRAME2 >= 0                                              (free superposition is a min)
    """
    from s15.phys_repl import averaged_backbone_from
    tg = I.targets()[:n_targets]
    a = b = c = 0.0
    fr_min = np.inf
    for t in tg:
        W, PHI, PSI, u = top75_windows(t["pdb"])
        W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
        nat = np.asarray(u["nat_ca"], float)
        for idx in (np.arange(len(W)), keep_lowest(np.arange(len(W))[::-1].astype(float), 37)):
            Y, Ybar, _C, _dev = common_frame_members(W[idx], PHI[idx], PSI[idx])
            avg, _C2, _d2 = averaged_backbone_from(W[idx], PHI[idx], PSI[idx])
            a = max(a, float(np.abs(Ybar - np.asarray(avg["CA"], float)).max()))
            k = kv_of(Y, nat, d_win=I.kabsch_rmsd_batch(W[idx], nat))
            b = max(b, abs(k["readout"] - I.ca_rmsd(Ybar, nat)))
            c = max(c, k["kv_resid"])
            fr_min = min(fr_min, k["FRAME2"])
    out = {"n_targets": int(n_targets),
           "max_abs_mean_delta_A": a, "max_readout_vs_metric_A": b,
           "max_kv_residual_A2": c, "min_FRAME2_A2": float(fr_min),
           "passed": bool(a < 1e-9 and b < 1e-9 and c < 1e-9 and fr_min > -1e-9)}
    if verbose:
        print(f"GC0  mean-delta {a:.3e} A | readout-vs-metric {b:.3e} A | "
              f"KV residual {c:.3e} A^2 | min FRAME^2 {fr_min:+.3e} A^2 -> "
              f"{'PASS' if out['passed'] else 'FAIL'}")
    return out


def gate_GC1(n_targets=8, verbose=True):
    """GC1 -- my Legacy scores ARE the Sprint-18 artefact's Legacy scores, candidate by candidate.

    `s18/results/down.json` stores `leg_terms` per candidate.  Recomputing them here through
    `s16.energy_lib.legacy_components_of_windows` must reproduce them exactly, or my gates are
    not the gates whose failure this sprint is explaining.
    """
    o = PL.read_complete(os.path.join(PL.RESULTS, "down.json"), need=126)
    rows = {r["pdb"]: r for r in o["rows"]}
    worst = 0.0
    worst_sp = 0.0
    for t in I.targets()[:n_targets]:
        r = rows[t["pdb"]]
        W, PHI, PSI, u = top75_windows(t["pdb"])
        comp, tot = legacy_scores(t["seq"], np.asarray(PHI, float), np.asarray(PSI, float))
        for k in LEG_TERMS:
            worst = max(worst, float(np.abs(comp[k] - np.asarray(r["leg_terms"][k], float)).max()))
        d = np.asarray(I.kabsch_rmsd_batch(np.asarray(W, float),
                                           np.asarray(u["nat_ca"], float)), float)
        worst_sp = max(worst_sp, float(np.abs(d - np.asarray(r["d"], float)).max()))
    out = {"n_targets": int(n_targets), "max_abs_leg_term_delta": worst,
           "max_abs_oracle_d_delta_A": worst_sp,
           "passed": bool(worst < 1e-9 and worst_sp < 1e-9)}
    if verbose:
        print(f"GC1  max |Legacy term delta| {worst:.3e} | max |ORACLE d delta| "
              f"{worst_sp:.3e} A -> {'PASS' if out['passed'] else 'FAIL'}")
    return out


if __name__ == "__main__":
    print("AGENT C gates -- these run before any scientific number is quoted.\n")
    g0 = gate_GC0()
    g1 = gate_GC1()
    PL.write("agentC_gates", {"GC0": g0, "GC1": g1,
                              "passed": bool(g0["passed"] and g1["passed"])},
             subdir="")
    import json
    with open(os.path.join(RESULTS, "agentC_gates.json"), "w") as fh:
        json.dump({"GC0": g0, "GC1": g1, "complete": True,
                   "passed": bool(g0["passed"] and g1["passed"])}, fh, indent=1)
    print(f"\nBOTH GATES {'PASS' if g0['passed'] and g1['passed'] else 'FAIL'}")
