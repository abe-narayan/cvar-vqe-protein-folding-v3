"""SPRINT 23, WORKSTREAM B -- H5: Ca-PRESERVING GEOMETRY REPAIR.

See `s23/PREREG_B.md` for the six-axis operator fork, hypothesis and falsifier, registered
BEFORE this script was run. Summary of the mechanism (full derivation in the prereg):

`core.amber.RESTRAINED_BACKBONE = ("N", "CA", "C")` is a MODULE-LEVEL GLOBAL, read once at
`AmberHamiltonian` construction time (`_index_topology`), not a per-call argument. The
deployed/incumbent restrained minimisation therefore pins the whole rigid peptide-plane frame
(N, Ca, C together), never Ca alone. This module builds a SECOND Hamiltonian per target with
`RESTRAINED_BACKBONE` monkey-patched to `("CA",)`, restrained-atom indices frozen at that
scope, and sweeps the SAME restraint-constant ladder on both scopes so the trade-off between
Ca displacement and physical validity can be read off directly and the CA-only scope compared
against the incumbent's own N+CA+C scope at matched k.

BASIS: point cloud throughout (top-75 coordinate average, rebuilt to a full backbone via
`s15.phys_repl.averaged_backbone_from` so AMBER has bonded atoms to act on). Never compared
to a rebuilt/projected chain's RMSD in this file.

Native Ca-RMSD is read for EVALUATION ONLY (ORACLE), never used to choose a k. Every other
axis (clash count, cis fraction, Ramachandran-outlier fraction, omega deviation, bond-length
deviation, bond-angle deviation, Ca displacement) is native-free.

AMBER/OpenMM SERIALISATION: this workstream is the sole authorised lane this sprint (BRIEF
Sec on AMBER). No live inter-agent channel was reachable from this session (no `ListAgents`),
so start/release is announced in the accompanying chat message and in `agentB_FINDINGS.md`.

Run:
    python -m s23.agentB_h5_carestraint --mode smoke     # 1 target, 2 k values, both scopes
    python -m s23.agentB_h5_carestraint --mode sweep      # full pre-registered n=30 sweep
    python -m s23.agentB_h5_carestraint --mode report     # summarise whatever has been written
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

from s12 import instrument as I                                   # noqa: E402
from s14.avgspace import top75_windows                            # noqa: E402
from s14.ener_avgrefine import geom_of, rama_ok, torsions_of, ATOMS  # noqa: E402
from s15.phys_repl import averaged_backbone_from                  # noqa: E402
from s16.repair import shape_subsample                            # noqa: E402
from core.amber import CONVERGE_MAX_KCAL, K_WEAK, K_MODERATE, K_STRONG  # noqa: E402

K_LADDER = (1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0)
SCOPES = ("ca", "full")
TOL = 1.0
STEPS = 0

IDEAL_BOND = dict(N_CA=1.458, CA_C=1.525, C_N=1.329)
IDEAL_ANGLE = dict(ang_N_CA_C=111.0, ang_CA_C_N=116.2, ang_C_N_CA=121.7)


# ------------------------------------------------------------------ atomic write
def _atomic_write(path, obj):
    d = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".tmp_", suffix=".json")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(obj, fh, indent=1,
                      default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


# ------------------------------------------------------------------ builders (bypass builder_for cache)
def _make_hamiltonian(seq, rep, scope):
    """A fresh AmberHamiltonian with RESTRAINED_BACKBONE set to `scope`'s atom set.

    Bypasses `core.amber.builder_for`/`_BUILDERS` deliberately: that cache keys on
    (sequence, n_bits, n_states, platform, threads) only, NOT on restraint scope, so reusing
    it across scopes would silently serve a builder frozen with the wrong `_restraint_idx`.
    """
    import core.amber as A
    old = A.RESTRAINED_BACKBONE
    A.RESTRAINED_BACKBONE = ("CA",) if scope == "ca" else ("N", "CA", "C")
    try:
        H = A.AmberHamiltonian(seq, rep, restraint_k=A.K_MODERATE, minimization_steps=0,
                               platform_name="CPU", collapse_floor=float("-inf"), threads=1)
    finally:
        A.RESTRAINED_BACKBONE = old
    assert scope != "ca" or H._restraint_idx and len(H._restraint_idx) == len(seq), (
        "CA-only restraint index count must equal residue count")
    assert scope != "full" or len(H._restraint_idx) == 3 * len(seq), (
        "N+CA+C restraint index count must equal 3x residue count")
    return H


def _drop(H):
    import core.amber as A
    A._drop_builder(H)


def _run_on(H, coords, k_restraint, components=True):
    import core.amber as A
    heavy = H._heavy_positions(coords, chi1=None)
    return A._run_memo(H, heavy, float(k_restraint), STEPS, TOL, components, memo=False)


# ------------------------------------------------------------------ validity metrics
def full_heavy_clash(heavy, heavy_names, min_sep=2, thresh=2.0):
    heavy = np.asarray(heavy, float)
    res = np.array([r for r, _ in heavy_names], int)
    D = np.linalg.norm(heavy[:, None, :] - heavy[None, :, :], axis=-1)
    D = np.where(np.abs(res[:, None] - res[None, :]) >= min_sep, D, 9e9)
    return int((D < thresh).sum() // 2), float(D.min())


def bonded_metrics(bb):
    from s14.ener_geom import _ang, _dih
    N, CA, Cc = (np.asarray(bb[k], float) for k in ("N", "CA", "C"))
    b = dict(N_CA=np.linalg.norm(CA - N, axis=1), CA_C=np.linalg.norm(Cc - CA, axis=1),
             C_N=np.linalg.norm(N[1:] - Cc[:-1], axis=1))
    a = dict(ang_N_CA_C=_ang(N, CA, Cc), ang_CA_C_N=_ang(CA[:-1], Cc[:-1], N[1:]),
             ang_C_N_CA=_ang(Cc[:-1], N[1:], CA[1:]))
    bond_dev = float(np.sqrt(np.mean([((b[k].mean() - IDEAL_BOND[k]) / IDEAL_BOND[k]) ** 2
                                       for k in IDEAL_BOND])))
    angle_dev = float(np.sqrt(np.mean([((a[k].mean() - IDEAL_ANGLE[k]) / IDEAL_ANGLE[k]) ** 2
                                        for k in IDEAL_ANGLE])))
    om = np.degrees(_dih(CA[:-1], Cc[:-1], N[1:], CA[1:]))
    omega_dev_mean = float(np.mean(np.abs(np.abs(om) - 180.0)))
    omega_dev_max = float(np.max(np.abs(np.abs(om) - 180.0))) if len(om) else float("nan")
    cis_frac = float(np.mean(np.abs(om) < 90.0)) if len(om) else float("nan")
    return dict(bond_len_rel_dev=bond_dev, bond_angle_rel_dev=angle_dev,
                omega_dev_mean=omega_dev_mean, omega_dev_max=omega_dev_max, cis_frac=cis_frac)


def score_result(out, input_ca, nat, seq):
    """Every axis, kept SEPARATE. `ca_rmsd_native` is the only ORACLE row."""
    ca = np.asarray(out["ca"], float)
    bb = {a: np.asarray(out["backbone"][a], float) for a in ("N", "CA", "C", "O")
          if a in out["backbone"]}
    heavy, heavy_names = out["heavy"], out["heavy_names"]
    n_clash, min_heavy = full_heavy_clash(heavy, heavy_names)
    phi, psi = torsions_of(bb)
    rama_frac_ok = rama_ok(phi, psi)
    bm = bonded_metrics(bb)
    ca_disp = float(np.sqrt(np.mean(np.sum((ca - np.asarray(input_ca, float)) ** 2, axis=1))))
    return {
        "ca_rmsd_native_ORACLE_EVAL_ONLY": float(I.ca_rmsd(ca, nat)),
        "ca_displacement_from_input": ca_disp,
        "n_clash_heavy_2A": n_clash,
        "min_heavy_dist": min_heavy,
        "rama_favoured_frac": float(rama_frac_ok),
        "rama_outlier_frac": float(1.0 - rama_frac_ok) if np.isfinite(rama_frac_ok) else None,
        **bm,
        "energy": float(out["energy"]), "energy_initial": float(out["energy_initial"]),
        "converged": bool(out["converged"]), "converge_reason": out["converge_reason"],
        "restraint_rmsd_of_restrained_set": float(out["restraint_rmsd"]),
        "wall": float(out["wall"]),
    }


# ------------------------------------------------------------------ one target
def run_target(t):
    import torsion_lib2 as tl2
    pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
    n = len(seq)
    W, PHI, PSI, u = top75_windows(pdb)
    nat = np.asarray(u["nat_ca"], float)
    avg, C_ca, dev = averaged_backbone_from(W, PHI, PSI)
    coords = {a: np.asarray(avg[a], float) for a in ATOMS if a in avg}

    rec = dict(pdb=pdb, fold=fold, n=n, seq=seq,
               pointcloud_dev_avg_vs_Cca=float(dev),
               incumbent_pointcloud_rmsd=float(I.ca_rmsd(C_ca, nat)))

    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)

    # k = 0: restraint-scope-invariant, run once.
    Hf = _make_hamiltonian(seq, rep, "full")
    out0 = _run_on(Hf, coords, 0.0, components=True)
    rec["k0_none"] = score_result(out0, C_ca, nat, seq)

    for scope in SCOPES:
        H = Hf if scope == "full" else _make_hamiltonian(seq, rep, "ca")
        for k in K_LADDER:
            out = _run_on(H, coords, k, components=True)
            rec[f"scope_{scope}_k{k:g}"] = score_result(out, C_ca, nat, seq)
        if scope == "ca":
            _drop(H)
    _drop(Hf)
    return rec


# ------------------------------------------------------------------ driver
def _path(tag="sweep"):
    return os.path.join(RESULTS, f"agentB_h5_carestraint_{tag}.json")


def run_sweep(targets=None, resume=True, tag="sweep"):
    tg = targets if targets is not None else [t for t in I.targets()
                                               if t["pdb"] in set(shape_subsample(30))]
    p = _path(tag)
    rows, done = [], set()
    if resume and os.path.exists(p):
        with open(p) as fh:
            rows = json.load(fh).get("per_target", [])
        done = {r["pdb"] for r in rows}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        try:
            rec = run_target(t)
        except Exception as exc:                                         # noqa: BLE001
            import core.amber as A
            A.clear_cache()
            rec = dict(pdb=t["pdb"], fold=int(t["fold"]), n=len(t["seq"]), seq=t["seq"],
                       error=str(exc)[:500])
            print(f"[{t['pdb']}] ERROR: {exc}", flush=True)
        rows.append(rec)
        if "error" not in rec:
            print(f"[{len(rows)}/{len(tg)}] {rec['pdb']} n={rec['n']} "
                  f"incumbent={rec['incumbent_pointcloud_rmsd']:.3f} "
                  f"k0={rec['k0_none']['ca_rmsd_native_ORACLE_EVAL_ONLY']:.3f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
        else:
            print(f"[{len(rows)}/{len(tg)}] {rec['pdb']} ({time.time() - t0:.0f}s)", flush=True)
        _write(rows, len(tg), tg, tag)
    _write(rows, len(tg), tg, tag)
    return rows


def _write(rows, n_expected, tg, tag="sweep"):
    target_ids = sorted(t["pdb"] for t in tg)
    obj = dict(what="H5: Ca-preserving geometry repair, k-ladder x restraint scope",
              k_ladder=list(K_LADDER), scopes=list(SCOPES), steps=STEPS, tolerance=TOL,
              converge_max_kcal=float(CONVERGE_MAX_KCAL),
              landmarks=dict(K_WEAK=K_WEAK, K_MODERATE=K_MODERATE, K_STRONG=K_STRONG),
              target_ids=target_ids, per_target=rows, n_rows=len(rows),
              n_expected=int(n_expected), complete=bool(len(rows) >= int(n_expected)))
    _atomic_write(_path(tag), obj)


# ------------------------------------------------------------------ paired stats (Hard Rules)
def stats(tag="sweep"):
    """Paired, fold-clustered-beside-iid comparisons on every axis, per BRIEF Hard Rules.

    Reuses `s23.qc_lib.paired_stats`/`verdict` -- Workstream C's shared, pre-registered
    fold-clustered-bootstrap implementation (`s23/PREREG_C.md`) -- rather than a second,
    divergent implementation of the same statistic.

    (1) scope_ca vs scope_full AT MATCHED k -- does restricting the restrained set change
        the trade-off, holding the restraint constant fixed (the H5 hypothesis directly).
    (2) each scope's own k=0 baseline -- the accuracy/validity cost of turning ANY restraint
        on, both scopes separately.
    """
    from s23 import qc_lib as QC
    p = _path(tag)
    if not os.path.exists(p):
        print("no results yet")
        return None
    with open(p) as fh:
        d = json.load(fh)
    rows = [r for r in d["per_target"] if "error" not in r]
    if len(rows) < 4:
        print(f"only {len(rows)} usable rows -- too few for a paired CI, skipping stats")
        return None
    pdbs = [r["pdb"] for r in rows]
    fold_of = {t["pdb"]: int(t["fold"]) for t in I.targets()}
    folds = np.array([fold_of[p_] for p_ in pdbs], int)

    axes = ("ca_rmsd_native_ORACLE_EVAL_ONLY", "ca_displacement_from_input",
            "n_clash_heavy_2A", "cis_frac", "rama_outlier_frac", "omega_dev_mean",
            "bond_len_rel_dev", "bond_angle_rel_dev")

    def col(key, ax):
        return np.array([r[key][ax] for r in rows], float)

    out = {"n": len(rows), "pdbs": pdbs, "comparisons": {}}
    for k in d["k_ladder"]:
        ka, kf = f"scope_ca_k{k:g}", f"scope_full_k{k:g}"
        if ka not in rows[0] or kf not in rows[0]:
            continue
        for ax in axes:
            a, b = col(ka, ax), col(kf, ax)
            label = f"ca_minus_full__k{k:g}__{ax}"
            c = QC.paired_stats(a, b, folds, label=label)
            c["verdict"] = QC.verdict(c)
            out["comparisons"][label] = c
    for scope in d["scopes"]:
        k0 = col("k0_none", "ca_rmsd_native_ORACLE_EVAL_ONLY")
        for k in d["k_ladder"]:
            key = f"scope_{scope}_k{k:g}"
            if key not in rows[0]:
                continue
            a = col(key, "ca_rmsd_native_ORACLE_EVAL_ONLY")
            label = f"{scope}_k{k:g}_minus_k0__ca_rmsd"
            c = QC.paired_stats(a, k0, folds, label=label)
            c["verdict"] = QC.verdict(c)
            out["comparisons"][label] = c
    _atomic_write(os.path.join(RESULTS, f"agentB_h5_carestraint_{tag}_stats.json"), out)
    for name, c in out["comparisons"].items():
        print(f"{name:<48} {c['mean_diff']:+.4f} SE {c['se']:.4f} MDE {c['mde']:.4f} "
              f"iid[{c['ci95_iid'][0]:+.4f},{c['ci95_iid'][1]:+.4f}] "
              f"fold[{c['ci95_fold'][0]:+.4f},{c['ci95_fold'][1]:+.4f}] "
              f"{c['n_better']}W/{c['n_worse']}L worst {c['worst_degradation']:+.3f}  "
              f"{c['verdict']}")
    return out


# ------------------------------------------------------------------ report
def report(tag="sweep"):
    p = _path(tag)
    if not os.path.exists(p):
        print("no results yet")
        return
    with open(p) as fh:
        d = json.load(fh)
    rows = [r for r in d["per_target"] if "error" not in r]
    errs = [r for r in d["per_target"] if "error" in r]
    print(f"complete={d['complete']} n_rows={d['n_rows']}/{d['n_expected']} "
          f"ok={len(rows)} errors={len(errs)}")
    for e in errs:
        print(f"  ERROR {e['pdb']}: {e.get('error')}")
    if not rows:
        return

    def col(key, sub):
        return np.array([r[key][sub] for r in rows], float)

    inc = np.array([r["incumbent_pointcloud_rmsd"] for r in rows], float)
    print(f"\nincumbent point cloud (no repair at all) Ca-RMSD: mean {inc.mean():.3f} n={len(inc)}")

    k0 = col("k0_none", "ca_rmsd_native_ORACLE_EVAL_ONLY")
    print(f"k=0 (unrestrained AMBER minimisation) Ca-RMSD: mean {k0.mean():.3f}  "
          f"clash {col('k0_none','n_clash_heavy_2A').mean():.2f}  "
          f"cis {col('k0_none','cis_frac').mean():.3f}  "
          f"rama_outlier {col('k0_none','rama_outlier_frac').mean():.3f}  "
          f"omega_dev {col('k0_none','omega_dev_mean').mean():.2f}  "
          f"bond_dev {col('k0_none','bond_len_rel_dev').mean():.4f}  "
          f"angle_dev {col('k0_none','bond_angle_rel_dev').mean():.4f}")

    print(f"\n{'scope':<6}{'k':>8}{'CaRMSD':>9}{'CaDisp':>9}{'clash':>8}{'cis':>7}"
          f"{'ramaOut':>9}{'omegaDv':>9}{'bondDv':>9}{'angDv':>9}{'conv%':>7}")
    K_LADDER_ = d["k_ladder"]
    for scope in d["scopes"]:
        for k in K_LADDER_:
            key = f"scope_{scope}_k{k:g}"
            if key not in rows[0]:
                continue
            rm = col(key, "ca_rmsd_native_ORACLE_EVAL_ONLY")
            dp = col(key, "ca_displacement_from_input")
            cl = col(key, "n_clash_heavy_2A")
            ci = col(key, "cis_frac")
            ro = col(key, "rama_outlier_frac")
            od = col(key, "omega_dev_mean")
            bd = col(key, "bond_len_rel_dev")
            ad = col(key, "bond_angle_rel_dev")
            cv = np.array([r[key]["converged"] for r in rows], bool)
            print(f"{scope:<6}{k:>8g}{rm.mean():>9.3f}{dp.mean():>9.3f}{cl.mean():>8.2f}"
                  f"{ci.mean():>7.3f}{ro.mean():>9.3f}{od.mean():>9.2f}{bd.mean():>9.4f}"
                  f"{ad.mean():>9.4f}{100*cv.mean():>6.0f}%")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="sweep", choices=("smoke", "sweep", "report", "stats"))
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "report":
        report()
    elif a.mode == "stats":
        stats()
    elif a.mode == "smoke":
        K_LADDER = (10.0, 100.0)
        tg = [t for t in I.targets() if t["pdb"] in set(shape_subsample(30))][:1]
        rows = run_sweep(targets=tg, resume=False, tag="smoke")
        print(json.dumps(rows, indent=1, default=lambda o: o.tolist()
                          if hasattr(o, "tolist") else str(o))[:4000])
    else:
        tg = [t for t in I.targets() if t["pdb"] in set(shape_subsample(30))]
        if a.n:
            tg = tg[:a.n]
        run_sweep(targets=tg, tag="sweep")
