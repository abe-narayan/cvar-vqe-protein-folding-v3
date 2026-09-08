"""s20/c_q1.py -- Q1: DO LEGACY AND AMBER CONTAIN COMPLEMENTARY INFORMATION?

Pre-registration: `s20/PREREG_C.md` section 2, written before this module produced a number.
Falsifier F-C1 registered there, first.

THE QUESTION, AND WHAT IT IS NOT.  Sprint 19 closed BOTH potentials in every decision role
(objective term, ranker, gate, steric rejector).  This module does not reopen any of them and
pre-commits that no result here is used to build a gate, a ranker or a rejector.  The open
question is what each model KNOWS -- which does not require either to be a good ranker.

THE INSTRUMENT.  126 targets x the SAME shipped top-75 windows.  On each candidate's
ideal-geometry rebuild:

    E_Legacy    genuine `core.energy`, eleven components, DEFAULT_WEIGHTS, never fitted
    E_AMBER     genuine ff14SB/GBn2 single point (`s17.phys_lib.ConstrainedBox.energy_point`)
    E_disto     the shipped leave-fold-out Bayes-risk distogram score
    d_reb       ORACLE Ca-RMSD of the REBUILD (not the window) to the native
    panel       `s16.energy_lib.panel`, kept as a VECTOR, never fused into a scalar
    align       ORACLE <e_i, b_pool>/||b_pool|| in the s19 common frame

Both energies are read from `s18/results/down.json` (complete, 126 rows), which computed them on
the identical rebuilds.  Gate GC20a reproduces the eleven Legacy components from scratch and
requires 0.00e+00 against that artefact before any Q1 number is quoted -- which certifies that
`s_amber_sp` is the AMBER single point of MY candidates and not of some other set.

BASIS (BRIEF section 1).  The ablation table is on the POINT-CLOUD basis: its readout is the
deployed coordinate average, comparable to Sprint 19's 3.050 A column and to nothing else.  Its
validity vector is reported on the MEMBERS (built chains), because a 22%-contracted point cloud
has no meaningful bond strain.

    python -m s20.c_q1 --gate       # GC20a
    python -m s20.c_q1 --smoke      # 6 targets, writes _SMOKE_*, never a result
    python -m s20.c_q1              # n = 126
    python -m s20.c_q1 --report
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                     # noqa: E402
from s14.avgspace import top75_windows              # noqa: E402
from s15 import seed as SD                          # noqa: E402
from s16 import energy_lib as EL                    # noqa: E402
from s18 import phys_lib as PL                      # noqa: E402
from s19 import agentC_lib as CL                    # noqa: E402

#: the fraction of the pool each model is said to "prefer".  1/3 -> 25 of 75, so the
#: independence expectation for the BOTH cell is 75/9 = 8.33 and the two disagreement cells are
#: 16.67 each.  Declared before the run; not swept.
PREF_Q = 1.0 / 3.0
M_KEEP = 38            # the Sprint-19 gate count, unchanged
N_RAND = 30            # matched-random draws for the ablation table
N_PART_NULL = 200      # matched-random partitions for the partition contrast null

#: the panel axes that are NOT constant by construction on an ideal-geometry rebuild.  bond
#: strain, angle strain, omega/cis and chirality are exactly ideal for every candidate (measured:
#: bond_strain 5.2e-06, angle_strain 2.4e-16, cis_frac 0.0, chirality_L_frac 1.0), so they carry
#: no information HERE and are reported as constant rather than as a null.
PANEL_LIVE = ("min_heavy", "n_clash_2A", "n_clash_2p6A",
              "rama_favoured", "rama_allowed", "rama_outlier")
PANEL_CONST = ("bond_strain", "angle_strain", "cis_frac", "chirality_L_frac")

AXES = ("d_reb", "rg", "align", "e_disto", "z_leg", "z_amb") + PANEL_LIVE


def _z(x):
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / s if s > 1e-12 else np.zeros_like(x)


def _spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 4:
        return float("nan")
    ra = np.argsort(np.argsort(a[ok])).astype(float)
    rb = np.argsort(np.argsort(b[ok])).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


def _rg(ca):
    ca = np.asarray(ca, float)
    return float(np.sqrt(((ca - ca.mean(0)) ** 2).sum(1).mean()))


# ==========================================================================
# per-target computation
# ==========================================================================
def load_down():
    o = PL.read_complete(os.path.join(PL.RESULTS, "down.json"), need=126)
    return {r["pdb"]: r for r in o["rows"]}


def target_row(t, drow, want_ablation=True):
    from core import geometry as geo
    pdb, seq, n, fold = t["pdb"], t["seq"], int(t["n"]), int(t["fold"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    nat = np.asarray(u["nat_ca"], float)
    K = len(W)

    #: the two genuine energies, on the identical rebuilds, read from the s18 artefact.
    comp = {k: np.asarray(v, float) for k, v in drow["leg_terms"].items()}
    e_leg = np.asarray(EL.legacy_total_from(comp), float)
    e_amb = np.asarray(drow["s_amber_sp"], float)

    #: the deployed operator's common frame -- members, mean, and the ORACLE error field.
    Y, Ybar, _C, _dev = CL.common_frame_members(W, PHI, PSI)
    T = I.superpose_batch(nat[None], Ybar)[0]
    e_field = Y - T                                          # (K, n, 3)  ORACLE
    b_pool = e_field.mean(0)
    nb = float(np.sqrt((b_pool ** 2).sum()))
    align = ((e_field.reshape(K, -1) @ b_pool.reshape(-1)) / nb) / np.sqrt(n) if nb > 1e-12 \
        else np.full(K, np.nan)                              # ORACLE, per-residue-rms units
    d_reb = np.asarray(I.kabsch_rmsd_batch(Y, nat), float)   # ORACLE

    #: structure of each candidate
    bb = geo.build_backbone_batch(PHI, PSI)
    pan = []
    for b in range(K):
        pan.append(EL.panel({a: np.asarray(bb[a][b], float) for a in bb}, seq))
    rg = np.array([_rg(bb["CA"][b]) for b in range(K)])

    #: the shipped distogram score on the same rebuilds (target-conditioned, native-free).
    dg = I.distogram(pdb, seq, fold)
    ii, jj = np.asarray(dg["i"], int), np.asarray(dg["j"], int)
    CAb = np.asarray(bb["CA"], float)
    D = np.linalg.norm(CAb[:, ii, :] - CAb[:, jj, :], axis=-1)
    e_dis = np.asarray(I.shipped_score(dg, D), float)

    ax = {"d_reb": d_reb, "rg": rg, "align": align, "e_disto": e_dis,
          "z_leg": _z(e_leg), "z_amb": _z(e_amb)}
    for k in PANEL_LIVE + PANEL_CONST:
        ax[k] = np.array([p[k] for p in pan], float)

    rec = {"pdb": pdb, "n": n, "fold": fold, "K": int(K),
           "rho_leg_amb": _spearman(e_leg, e_amb),
           "rho_leg_d": _spearman(e_leg, d_reb),
           "rho_amb_d": _spearman(e_amb, d_reb),
           "rho_leg_dis": _spearman(e_leg, e_dis),
           "rho_amb_dis": _spearman(e_amb, e_dis),
           "rho_dis_d": _spearman(e_dis, d_reb),
           "b_pool_norm_per_res": float(nb / np.sqrt(n)),
           "readout": float(I.ca_rmsd(Ybar, nat))}

    #: ---- THE PARTITION.  Cells are defined by RANK, so they are native-free.
    nq = max(1, int(round(K * PREF_Q)))
    top_l = set(CL.keep_lowest(e_leg, nq).tolist())
    top_a = set(CL.keep_lowest(e_amb, nq).tolist())
    cells = {"both": sorted(top_l & top_a),
             "leg_only": sorted(top_l - top_a),
             "amb_only": sorted(top_a - top_l),
             "neither": sorted(set(range(K)) - top_l - top_a)}
    rec["cell_sizes"] = {k: len(v) for k, v in cells.items()}
    rec["cells"] = {k: [int(x) for x in v] for k, v in cells.items()}
    rec["cell_mean"] = {c: {a: (float(np.nanmean(ax[a][idx])) if idx else float("nan"))
                            for a in AXES} for c, idx in cells.items()}

    #: matched-random partitions at IDENTICAL cell sizes -> the null for every contrast.
    rng = SD.stable_rng(pdb, "s20C_q1_part")
    sizes = [len(cells[c]) for c in ("both", "leg_only", "amb_only", "neither")]
    null = {a: [] for a in AXES}
    for _ in range(N_PART_NULL):
        p = rng.permutation(K)
        o = 0; parts = []
        for s in sizes:
            parts.append(p[o:o + s]); o += s
        for a in AXES:
            v = ax[a]
            null[a].append(float(np.nanmean(v[parts[1]])) - float(np.nanmean(v[parts[2]]))
                           if sizes[1] and sizes[2] else np.nan)
    rec["contrast_leg_minus_amb"] = {a: (rec["cell_mean"]["leg_only"][a]
                                         - rec["cell_mean"]["amb_only"][a]) for a in AXES}
    rec["contrast_null"] = {a: {"mean": float(np.nanmean(null[a])),
                                "sd": float(np.nanstd(null[a])),
                                "lo": float(np.nanpercentile(null[a], 2.5)),
                                "hi": float(np.nanpercentile(null[a], 97.5))} for a in AXES}
    rec["contrast_z"] = {a: (float((rec["contrast_leg_minus_amb"][a]
                                    - rec["contrast_null"][a]["mean"])
                                   / rec["contrast_null"][a]["sd"]))
                         if rec["contrast_null"][a]["sd"] > 1e-12 else float("nan")
                         for a in AXES}

    #: ---- THE MANDATED ABLATION TABLE (directive section 23).  POINT-CLOUD basis.
    if want_ablation:
        rec["ablation"] = ablation(W, PHI, PSI, nat, seq, e_leg, e_amb, e_dis, ax, pdb)
    return rec


def _emit(W, PHI, PSI, nat, idx):
    Y, Ybar, _C, _d = CL.common_frame_members(W[idx], PHI[idx], PSI[idx])
    D2 = float(np.sum((Y - Ybar) ** 2) / (Y.shape[0] * Y.shape[1]))
    return {"rmsd": float(I.ca_rmsd(Ybar, nat)), "D": float(np.sqrt(max(D2, 0.0)))}


def _arm(W, PHI, PSI, nat, idx, ax):
    o = _emit(W, PHI, PSI, nat, idx)
    o["best_member"] = float(np.min(ax["d_reb"][idx]))
    o["mean_member"] = float(np.mean(ax["d_reb"][idx]))
    for a in PANEL_LIVE + PANEL_CONST:
        o["v_" + a] = float(np.nanmean(ax[a][idx]))
    o["mean_align"] = float(np.nanmean(ax["align"][idx]))
    o["m"] = int(len(idx))
    return o


def ablation(W, PHI, PSI, nat, seq, e_leg, e_amb, e_dis, ax, pdb):
    """Random / Distogram / Legacy / AMBER / Legacy->AMBER / AMBER->Legacy, identical candidates,
    identical count, and a MATCHED-RANDOM ordering for the one-stage and two-stage arms alike."""
    K = len(W)
    m = M_KEEP
    m1 = int(round(K * np.sqrt(m / K)))     # 53 for K=75, m=38: two equal-ratio stages
    out = {}
    out["none"] = _arm(W, PHI, PSI, nat, np.arange(K), ax)
    for nm, s in (("legacy", e_leg), ("amber", e_amb), ("disto", e_dis)):
        out[nm] = _arm(W, PHI, PSI, nat, CL.keep_lowest(s, m), ax)
    for nm, (s1, s2) in (("leg_then_amb", (e_leg, e_amb)), ("amb_then_leg", (e_amb, e_leg))):
        i1 = CL.keep_lowest(s1, m1)
        i2 = i1[CL.keep_lowest(np.asarray(s2)[i1], m)]
        out[nm] = _arm(W, PHI, PSI, nat, np.sort(i2), ax)
    #: MATCHED-RANDOM, one stage and two stages.  Same count, same number of stages.
    rng = SD.stable_rng(pdb, "s20C_q1_abl")
    for nm, two in (("rand", False), ("rand2", True)):
        acc = []
        for _ in range(N_RAND):
            if two:
                i1 = np.sort(rng.permutation(K)[:m1])
                idx = np.sort(i1[rng.permutation(len(i1))[:m]])
            else:
                idx = np.sort(rng.permutation(K)[:m])
            acc.append(_arm(W, PHI, PSI, nat, idx, ax))
        out[nm] = {k: float(np.mean([a[k] for a in acc])) for k in acc[0]}
        out[nm]["n_draws"] = N_RAND
    return out


# ==========================================================================
# gate GC20a
# ==========================================================================
def gate_GC20a(n_targets=8, verbose=True):
    """My Legacy components ARE `s18/results/down.json`'s, candidate by candidate, which
    certifies that its `s_amber_sp` is the AMBER single point of MY candidates."""
    down = load_down()
    worst_c = worst_d = 0.0
    rows = []
    for t in I.targets()[:n_targets]:
        pdb, seq = t["pdb"], t["seq"]
        W, PHI, PSI, u = top75_windows(pdb)
        comp = EL.legacy_components_of_windows(seq, np.asarray(PHI, float),
                                               np.asarray(PSI, float))
        ref = down[pdb]["leg_terms"]
        dc = max(float(np.abs(np.asarray(comp[k], float)
                              - np.asarray(ref[k], float)).max()) for k in comp)
        dd = float(np.abs(np.asarray(I.kabsch_rmsd_batch(np.asarray(W, float),
                                                         np.asarray(u["nat_ca"], float)), float)
                          - np.asarray(down[pdb]["d"], float)).max())
        worst_c = max(worst_c, dc); worst_d = max(worst_d, dd)
        rows.append({"pdb": pdb, "d_terms": dc, "d_oracle": dd})
    out = {"n_targets": int(n_targets), "rows": rows,
           "max_d_terms": worst_c, "max_d_oracle_d": worst_d,
           "passed": bool(worst_c < 1e-9 and worst_d < 1e-9)}
    if verbose:
        print(f"GC20a  max |delta Legacy term| {worst_c:.3e} | max |delta ORACLE d| "
              f"{worst_d:.3e} A over {n_targets} targets -> "
              f"{'PASS' if out['passed'] else 'FAIL'}")
    return out


# ==========================================================================
def run(limit=0, out="c_q1.json", verbose=True):
    down = load_down()
    tg = I.targets()
    if limit:
        tg = tg[:limit]
    path = os.path.join(RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            rows = json.load(open(path)).get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    cfg = {"PREF_Q": PREF_Q, "M_KEEP": M_KEEP, "N_RAND": N_RAND,
           "N_PART_NULL": N_PART_NULL, "AXES": list(AXES),
           "source": "s18/results/down.json"}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        rows.append(target_row(t, down[t["pdb"]]))
        if verbose and len(rows) % 10 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
        _write(out, rows, cfg, len(tg))
    _write(out, rows, cfg, len(tg))
    print(f"DONE {len(rows)}/{len(tg)} in {time.time()-t0:.0f}s", flush=True)
    return rows


def _write(out, rows, cfg, n_expected):
    obj = {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg),
           "n_rows": len(rows), "n_expected": int(n_expected),
           "complete": bool(len(rows) >= int(n_expected))}
    p = os.path.join(RESULTS, out)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, p)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", action="store_true")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    if a.gate:
        g = gate_GC20a()
        with open(os.path.join(RESULTS, "c_q1_gate.json"), "w") as fh:
            json.dump({"GC20a": g, "complete": True, "passed": g["passed"]}, fh, indent=1)
        sys.exit(0 if g["passed"] else 1)
    if a.smoke:
        run(limit=6, out="_SMOKE_c_q1.json")
    else:
        run(limit=a.limit)
