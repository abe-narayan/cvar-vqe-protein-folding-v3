"""s22/probe_perturb.py -- WORKSTREAM B, PRIORITY 1: THE CONTROLLED SYNTHETIC PERTURBATION PROBE.

`s22/PREREG_B.md` sections 1-9, written and frozen BEFORE this file produced any number.

WHAT THIS MEASURES.  Five matched-magnitude torsion-space perturbation constructions ("axes") are
applied to real starting structures (the shipped top-75 real-rebuild pool, medoid + 2 pool draws per
target), and dE_Legacy / dE_AMBER are read at every trial.  AMBER here is the BARE SINGLE POINT
(`ConstrainedBox.energy_point` via `s20.c_land.Pot.amber`) -- see PREREG_B.md section 9, FUNCTIONAL
fork: a landscape mechanism probe needs a function of theta, and the deployed `H_AMBER = E o
Relax_50` is not one.  Every number in this file is about the BARE single point; it is restated in
every printed table so no reader has to remember it.

THE FIVE AXES, each a torsion-space direction of CONTROLLED MAGNITUDE (matched RMS over the full
active coordinate set, so "magnitude" means the same thing on every axis):

    A  diffuse/random   -- iid direction spread over every active (phi,psi)
    B  compactness       -- coherent GLOBAL bias toward the extended-strand reference
                            (-120, 130 deg, the same basin `core.energy._RAMA_BASINS[1]` names, and
                            the same reference `core.amber.AmberHamiltonian._reference_backbone` uses)
                            or toward the helical reference (-63, -42), sign = which
    C  contact density   -- the SAME coherent bias, restricted to the middle third of the chain only
    D  steric overlap    -- a numerically-computed direction that most changes the distance of the
                            single closest non-bonded (|i-j| >= 3) CB-CB pair in the START structure
    E  local geometry     -- the SAME per-coordinate random construction as axis A, but concentrated
                            entirely onto one residue's (phi, psi) instead of spread over the chain

Falsifiers F-B1 (concentration ratio) and F-B2 (compactness specialisation) are in PREREG_B.md
section 4 and are evaluated by `analyse()` below, unedited from what was registered.

    python -m s22.probe_perturb run        # the 30-target subset
    python -m s22.probe_perturb --smoke    # 2 targets
    python -m s22.probe_perturb analyse
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
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

from s12 import instrument as I                       # noqa: E402
from s14.avgspace import top75_windows                # noqa: E402
from s15 import seed as SD                             # noqa: E402
from s20 import c_land as CL                            # noqa: E402
from s21 import c_norm as CN                            # noqa: E402

SUBSET = CN.SUBSET                       # the SAME 30 targets that carry the C5 curvature numbers
N_STARTS = 3
MAGS = (0.02, 0.05, 0.10)                # radians, RMS over the full active coordinate set
N_REPEAT_RANDOM = 4                      # axes A, E
SIGNS = (1.0, -1.0)                      # axes B, C, D

#: the two reference basins, same numeric values as `core.energy._RAMA_BASINS`, restated here so
#: this module has no import-time dependency on a private name.
REF_HELIX = (np.radians(-63.0), np.radians(-42.0))
REF_EXTENDED = (np.radians(-120.0), np.radians(130.0))


# ==========================================================================================
def wrap(x):
    return (np.asarray(x, float) + np.pi) % (2.0 * np.pi) - np.pi


def rms(v):
    v = np.asarray(v, float)
    return float(np.sqrt(np.mean(v * v))) if v.size else 0.0


def unit_rms(v, n_active):
    v = np.asarray(v, float)
    r = rms(v)
    if r < 1e-15:
        return v
    return v / r


def concentration_kappa(delta_active):
    """Participation-ratio-style concentration index of a perturbation vector, folded to
    per-RESIDUE (combining each residue's (phi,psi) pair into one magnitude) so a kick on one
    residue's two coordinates reads as concentration on ONE unit, not two.  1/n_res = fully
    diffuse, 1.0 = a single residue carries the whole move."""
    n_res = len(delta_active) // 2 + 1     # active excludes phi[0] only; approx, see caller
    a = np.abs(delta_active)
    s1, s2 = a.sum(), (a * a).sum()
    if s2 <= 0:
        return float("nan")
    p = (s1 ** 2) / (max(len(delta_active), 1) * s2)
    return float(p)


# ==========================================================================================
class TargetGeom:
    """Per-target geometric proxies, independent of either energy model."""

    def __init__(self, seq):
        from core import energy as et
        from core import geometry as geo
        self.et = et
        self.geo = geo
        self.seq = seq
        self.n = len(seq)
        di, dj, sep = et.pair_index(self.n)
        self.m3 = sep >= 3
        self.di, self.dj = di[self.m3], dj[self.m3]

    def coords(self, phi, psi):
        return self.geo.build_backbone(phi, psi)

    def rg(self, CA):
        return float(self.geo.radius_of_gyration(np.asarray(CA, float)))

    def contact_frac(self, coords):
        CB = np.asarray(coords["CB"], float)
        d = np.linalg.norm(CB[self.di] - CB[self.dj], axis=1)
        return float((d < 8.5).mean()) if len(d) else float("nan")

    def steric(self, coords):
        return float(self.et.steric_term(coords, self.seq, rings=None))

    def cb_pair_dist(self, coords, i, j):
        CB = np.asarray(coords["CB"], float)
        return float(np.linalg.norm(CB[i] - CB[j]))

    def closest_pair(self, coords):
        CB = np.asarray(coords["CB"], float)
        d = np.linalg.norm(CB[self.di] - CB[self.dj], axis=1)
        k = int(np.argmin(d))
        return int(self.di[k]), int(self.dj[k]), float(d[k])


# ==========================================================================================
def build_directions(theta0, n, TG, active):
    """All FIVE directed constructions at unit RMS over the FULL active set (len(active)).

    Returns a dict axis -> callable(rng_or_none) -> unit-RMS active-length vector, except D which
    is precomputed once (deterministic, native-free -- uses only the start's own geometry).
    """
    n_active = len(active)
    phi0, psi0 = theta0[:n], theta0[n:]

    def ref_vec(pair):
        p, s = pair
        return np.concatenate([np.full(n, p), np.full(n, s)])

    dir_ext = wrap(ref_vec(REF_EXTENDED) - theta0)[active]
    dir_hel = wrap(ref_vec(REF_HELIX) - theta0)[active]
    dir_ext_u = unit_rms(dir_ext, n_active)
    dir_hel_u = unit_rms(dir_hel, n_active)

    #: axis C -- restrict the SAME coherent bias to the middle third of the chain (phi AND psi of
    #: those residues), zero elsewhere, THEN renormalise to unit RMS over the FULL active length --
    #: this is what concentrates the identical torsion budget onto fewer residues.
    lo, hi = n // 3, n - n // 3
    window = set(range(lo, hi)) if hi > lo else {n // 2}
    local_mask_full = np.zeros(2 * n, dtype=bool)
    for r in window:
        local_mask_full[r] = True          # phi_r
        local_mask_full[n + r] = True       # psi_r
    local_mask_active = local_mask_full[active]

    def local_variant(dir_full_active):
        v = np.where(local_mask_active, dir_full_active, 0.0)
        return unit_rms(v, n_active)

    dir_ext_local = local_variant(dir_ext)
    dir_hel_local = local_variant(dir_hel)

    #: axis D -- numeric gradient of the closest non-bonded CB-CB distance w.r.t. every active
    #: coordinate, by central difference on the PURE GEOMETRY (no energy call).  Deterministic,
    #: native-free: the pair is chosen from the start structure's own closeness ranking.
    coords0 = TG.coords(phi0, psi0)
    i_pair, j_pair, d0_pair = TG.closest_pair(coords0)
    h = 1e-4
    grad_d = np.zeros(n_active)
    for a_, k in enumerate(active):
        tp = theta0.copy(); tp[k] += h
        tm = theta0.copy(); tm[k] -= h
        cp = TG.coords(tp[:n], tp[n:])
        cm = TG.coords(tm[:n], tm[n:])
        dp = TG.cb_pair_dist(cp, i_pair, j_pair)
        dm = TG.cb_pair_dist(cm, i_pair, j_pair)
        grad_d[a_] = (dp - dm) / (2.0 * h)
    #: sign convention: POSITIVE sign of the axis DECREASES d_ij (drives the pair together, "clash"
    #: direction); NEGATIVE sign increases it ("relief").  -grad_d moves d_ij down fastest.
    dir_steric_u = unit_rms(-grad_d, n_active)

    return {
        "ext_u": dir_ext_u, "hel_u": dir_hel_u,
        "ext_local": dir_ext_local, "hel_local": dir_hel_local,
        "steric_u": dir_steric_u,
        "pair": (i_pair, j_pair, d0_pair),
        "window": sorted(window),
    }


def single_residue_dir(n, active, mid, rng):
    """Axis E: a random unit-RMS direction over the FULL active set, but with all mass on
    residue `mid`'s (phi, psi) pair -- the concentration extreme opposite of axis A."""
    n_active = len(active)
    v = np.zeros(n_active)
    idx_phi = int(np.where(active == mid)[0][0]) if mid in active else None
    idx_psi = int(np.where(active == (n + mid))[0][0])
    g = rng.normal(size=2)
    g = g / max(np.linalg.norm(g), 1e-12)
    if idx_phi is not None:
        v[idx_phi] = g[0]
    v[idx_psi] = g[1]
    return unit_rms(v, n_active)


def random_dir(n_active, rng):
    g = rng.normal(size=n_active)
    return unit_rms(g, n_active)


# ==========================================================================================
def eval_trial(P, TG, theta0, coords0, e_l0, e_a0, rg0, cf0, st0, seq, n, active,
               axis, sign, mag, repeat, direction, pair=None):
    theta1 = theta0.copy()
    theta1[active] = theta0[active] + mag * sign * direction
    theta1 = wrap(theta1)
    phi1, psi1 = theta1[:n], theta1[n:]
    coords1 = TG.coords(phi1, psi1)
    e_l1 = float(P.legacy(theta1)[0])
    e_a1 = float(P.amber(theta1)[0])
    rg1 = TG.rg(coords1["CA"])
    cf1 = TG.contact_frac(coords1)
    st1 = TG.steric(coords1)
    realised_rms = rms(wrap(theta1 - theta0)[active])
    kappa = concentration_kappa(mag * sign * direction)
    d_pair0 = d_pair1 = float("nan")
    if pair is not None:
        i_p, j_p, _ = pair
        d_pair0 = TG.cb_pair_dist(coords0, i_p, j_p)
        d_pair1 = TG.cb_pair_dist(coords1, i_p, j_p)
    return {
        "axis": axis, "sign": float(sign), "mag": float(mag), "repeat": int(repeat),
        "rms_realised": realised_rms, "kappa": kappa,
        "dRg": rg1 - rg0, "dcontact": cf1 - cf0, "dsteric": st1 - st0,
        "d_pair": (d_pair1 - d_pair0) if pair is not None else float("nan"),
        "dE_legacy": e_l1 - e_l0, "dE_amber": e_a1 - e_a0,
    }


def target_row(t, verbose=True):
    pdb, seq, n = t["pdb"], t["seq"], int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    W = np.asarray(W, float); PHI = np.asarray(PHI, float); PSI = np.asarray(PSI, float)
    Pm = I.pairwise_rmsd(W)
    med = int(I.medoid(Pm))
    rng0 = SD.stable_rng(pdb, "s20C_land")
    others = [int(x) for x in rng0.permutation(len(W))[:N_STARTS] if int(x) != med]
    order = [med] + others[:N_STARTS - 1]
    starts = [np.concatenate([PHI[b], PSI[b]]) for b in order]

    P = CL.Pot(seq)
    TG = TargetGeom(seq)
    active = np.arange(1, 2 * n)
    mid = max(1, n // 2)
    rec = {"pdb": pdb, "n": n, "fold": int(t["fold"]), "start_order": order, "trials": []}

    for si, theta0 in enumerate(starts):
        rngA = SD.stable_rng(pdb, f"s22probeA_{si}")
        rngE = SD.stable_rng(pdb, f"s22probeE_{si}")
        coords0 = TG.coords(theta0[:n], theta0[n:])
        e_l0 = float(P.legacy(theta0)[0])
        e_a0 = float(P.amber(theta0)[0])
        rg0 = TG.rg(coords0["CA"])
        cf0 = TG.contact_frac(coords0)
        st0 = TG.steric(coords0)
        dirs = build_directions(theta0, n, TG, active)

        common = dict(P=P, TG=TG, theta0=theta0, coords0=coords0, e_l0=e_l0, e_a0=e_a0,
                      rg0=rg0, cf0=cf0, st0=st0, seq=seq, n=n, active=active)

        # axis A -- diffuse random
        for m in MAGS:
            for r in range(N_REPEAT_RANDOM):
                d = random_dir(len(active), rngA)
                rec["trials"].append(eval_trial(**common, axis="A_random", sign=1.0, mag=m,
                                                 repeat=r, direction=d))
        # axis E -- single-residue concentrated
        for m in MAGS:
            for r in range(N_REPEAT_RANDOM):
                d = single_residue_dir(n, active, mid, rngE)
                rec["trials"].append(eval_trial(**common, axis="E_local_residue", sign=1.0, mag=m,
                                                 repeat=r, direction=d))
        # axis B -- compactness (global coherent bias); sign +1 = toward extended, -1 = toward helix
        for m in MAGS:
            for s in SIGNS:
                d = dirs["ext_u"] if s > 0 else dirs["hel_u"]
                rec["trials"].append(eval_trial(**common, axis="B_compactness", sign=s, mag=m,
                                                 repeat=0, direction=d))
        # axis C -- contact density (same bias, local window only)
        for m in MAGS:
            for s in SIGNS:
                d = dirs["ext_local"] if s > 0 else dirs["hel_local"]
                rec["trials"].append(eval_trial(**common, axis="C_contact_local", sign=s, mag=m,
                                                 repeat=0, direction=d))
        # axis D -- steric pair-directed; sign +1 = clash (drive together), -1 = relief
        for m in MAGS:
            for s in SIGNS:
                rec["trials"].append(eval_trial(**common, axis="D_steric_pair", sign=s, mag=m,
                                                 repeat=0, direction=dirs["steric_u"],
                                                 pair=dirs["pair"]))
        rec.setdefault("meta", []).append({
            "start_idx": si, "pair": dirs["pair"], "window": dirs["window"],
            "mid_residue": mid, "E_legacy_0": e_l0, "E_amber_0": e_a0, "Rg_0": rg0,
        })
        if verbose:
            print(f"    start {si}: {len(rec['trials'])} trials so far", flush=True)
    P.close()
    return rec


# ==========================================================================================
def run(subset=None, out=None, verbose=True):
    sub = SUBSET if subset is None else subset
    tg = {t["pdb"]: t for t in I.targets()}
    cfg = {"MAGS": list(MAGS), "N_STARTS": N_STARTS, "N_REPEAT_RANDOM": N_REPEAT_RANDOM,
           "SIGNS": list(SIGNS), "REF_HELIX": list(REF_HELIX), "REF_EXTENDED": list(REF_EXTENDED),
           "subset": list(sub), "amber_object": "bare single point (ConstrainedBox.energy_point), "
           "NO minimisation -- see PREREG_B.md section 9, FUNCTIONAL fork"}
    cfg_hash = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:16]
    out = out or f"probe_perturb_{cfg_hash}.json"
    path = os.path.join(RESULTS, out)
    rows, done = [], set()
    if os.path.exists(path):
        try:
            prev = json.load(open(path))
            if prev.get("cfg_hash") == cfg_hash:
                rows = prev.get("rows", [])
                done = {r["pdb"] for r in rows}
        except Exception:
            rows = []
    t0 = time.time()
    for pdb in sub:
        if pdb in done:
            continue
        if verbose:
            print(f"[{len(rows)+1}/{len(sub)}] {pdb} ...", flush=True)
        rows.append(CN.with_mem_retry(target_row, tg[pdb], verbose=verbose))
        _write(path, rows, cfg, cfg_hash, sub)
        if verbose:
            print(f"  done {pdb} ({time.time()-t0:.0f}s elapsed)", flush=True)
    _write(path, rows, cfg, cfg_hash, sub)
    print(f"DONE {len(rows)}/{len(sub)} in {time.time()-t0:.0f}s -> {path}", flush=True)
    return path


def _write(path, rows, cfg, cfg_hash, sub):
    #: completion demands the FULL key set per trial, not merely a row count (BRIEF section 3.5).
    required_trial_keys = {"axis", "sign", "mag", "repeat", "rms_realised", "kappa", "dRg",
                            "dcontact", "dsteric", "d_pair", "dE_legacy", "dE_amber"}
    full = True
    for r in rows:
        if not r.get("trials"):
            full = False; break
        for tr in r["trials"]:
            if not required_trial_keys.issubset(tr.keys()):
                full = False; break
        if not full:
            break
    obj = {"rows": rows, "config": cfg, "cfg_hash": cfg_hash,
           "n_rows": len(rows), "n_expected": len(sub),
           "complete": bool(len(rows) >= len(sub) and full and list({r['pdb'] for r in rows})
                            and set(r["pdb"] for r in rows) == set(sub))}
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    os.replace(tmp, path)


# ==========================================================================================
# ANALYSIS -- F-B1 (concentration ratio) and F-B2 (compactness specialisation)
# ==========================================================================================
AXES = ("A_random", "B_compactness", "C_contact_local", "D_steric_pair", "E_local_residue")
PROXY_OF = {"A_random": "rms_realised", "B_compactness": "dRg", "C_contact_local": "dcontact",
            "D_steric_pair": "d_pair", "E_local_residue": "rms_realised"}


def _load_norm():
    p = os.path.join(ROOT, "s21", "results", "c_norm.json")
    o = json.load(open(p))
    return {r["pdb"]: r for r in o["rows"]}


def _flatten(rows, normrow):
    """One row per trial, tagged with its target and the robust z of dE using the ALREADY-AUDITED
    per-target median/MAD from s21/results/c_norm.json (bare AMBER single point, genuine Legacy
    total, over that target's own top-75 pool) -- PREREG_B.md section 9, NORMALISATION fork."""
    out = []
    for r in rows:
        pdb = r["pdb"]
        nr = normrow.get(pdb)
        if nr is None:
            continue
        mL, sL = nr["legacy"]["median"], max(nr["legacy"]["mad_n"], 1e-30)
        mA, sA = nr["amber"]["median"], max(nr["amber"]["mad_n"], 1e-30)
        for tr in r["trials"]:
            row = dict(tr)
            row["pdb"] = pdb
            row["zL"] = tr["dE_legacy"] / sL
            row["zA"] = tr["dE_amber"] / sA
            out.append(row)
    return out


def _spearman(x, y):
    from scipy.stats import spearmanr
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5 or np.std(x[m]) < 1e-15 or np.std(y[m]) < 1e-15:
        return float("nan")
    return float(spearmanr(x[m], y[m]).correlation)


def _bootstrap_over_targets(trials, fn, n_boot=2000, seed=0):
    pdbs = sorted(set(t["pdb"] for t in trials))
    by_pdb = {p: [t for t in trials if t["pdb"] == p] for p in pdbs}
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(pdbs, size=len(pdbs), replace=True)
        sample = [t for p in pick for t in by_pdb[p]]
        v = fn(sample)
        if np.isfinite(v):
            vals.append(v)
    return np.asarray(vals, float)


def _sensitivity(trials, axis):
    sub = [t for t in trials if t["axis"] == axis]
    proxy = PROXY_OF[axis]
    x = [t[proxy] for t in sub]
    return {
        "n": len(sub),
        "spearman_legacy": _spearman(x, [t["zL"] for t in sub]),
        "spearman_amber": _spearman(x, [t["zA"] for t in sub]),
    }


def _rate(trials, axes):
    """median(|dE| / rms_realised) over trials whose axis is in `axes`, separately per potential.
    A ratio WITHIN one potential (units cancel), never compared raw across potentials."""
    sub = [t for t in trials if t["axis"] in axes and t["rms_realised"] > 1e-9]
    if not sub:
        return float("nan"), float("nan")
    rL = np.median([abs(t["dE_legacy"]) / t["rms_realised"] for t in sub])
    rA = np.median([abs(t["dE_amber"]) / t["rms_realised"] for t in sub])
    return float(rL), float(rA)


def analyse(path=None):
    if path is None:
        cands = sorted(f for f in os.listdir(RESULTS) if f.startswith("probe_perturb_"))
        if not cands:
            raise SystemExit("no probe_perturb_*.json in s22/results -- run first")
        path = os.path.join(RESULTS, cands[-1])
    o = json.load(open(path))
    if not o.get("complete"):
        print(f"*** WARNING: {path} is NOT complete ({o['n_rows']}/{o['n_expected']}) ***")
    normrow = _load_norm()
    trials = _flatten(o["rows"], normrow)

    print("=" * 100)
    print(f"PROBE ANALYSIS  n_targets={len(o['rows'])}  n_trials={len(trials)}  "
          f"source={os.path.basename(path)}")
    print("  AMBER = BARE SINGLE POINT throughout (ConstrainedBox.energy_point, NO minimisation).")
    print("=" * 100)

    sens = {ax: _sensitivity(trials, ax) for ax in AXES}
    print(f"\n{'axis':<20}{'proxy':<14}{'n':>6}{'spearman(Legacy)':>20}{'spearman(AMBER)':>18}")
    for ax in AXES:
        s = sens[ax]
        print(f"{ax:<20}{PROXY_OF[ax]:<14}{s['n']:>6}{s['spearman_legacy']:>20.4f}"
              f"{s['spearman_amber']:>18.4f}")

    # ---- F-B1: concentration ratio ----
    rL_conc, rA_conc = _rate(trials, ("D_steric_pair", "E_local_residue"))
    rL_diff, rA_diff = _rate(trials, ("A_random",))
    R_legacy = rL_conc / max(rL_diff, 1e-30)
    R_amber = rA_conc / max(rA_diff, 1e-30)
    ratio = R_amber / max(R_legacy, 1e-30)

    def _R_ratio_stat(sample):
        rl_c, ra_c = _rate(sample, ("D_steric_pair", "E_local_residue"))
        rl_d, ra_d = _rate(sample, ("A_random",))
        if rl_d <= 0 or rl_c != rl_c or ra_d <= 0:
            return float("nan")
        Rl = rl_c / rl_d
        Ra = ra_c / ra_d
        if Rl <= 0:
            return float("nan")
        return Ra / Rl

    boot = _bootstrap_over_targets(trials, _R_ratio_stat)
    ci = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))] if len(boot) else [
        float("nan")] * 2
    print(f"\nF-B1 CONCENTRATION RATIO  (median |dE|/rms over concentrated axes {{D,E}} "
          f"/ diffuse axis {{A}})")
    print(f"  R(Legacy) = {R_legacy:.4f}   R(AMBER) = {R_amber:.4f}")
    print(f"  R(AMBER)/R(Legacy) = {ratio:.4f}   bootstrap 95% CI (over {len(o['rows'])} targets, "
          f"{len(boot)} draws) = [{ci[0]:.4f}, {ci[1]:.4f}]")
    fb1_fires = bool(ci[0] > 1.0)
    fb1_refuted = bool(ci[1] < 1.0)
    print(f"  F-B1: {'SUPPORTED (CI > 1)' if fb1_fires else 'REFUTED (CI < 1, wrong side)' if fb1_refuted else 'NOT MEASURED (CI spans 1)'}")

    # ---- F-B2: compactness specialisation on axis B ----
    def _b2_stat(sample):
        s = _sensitivity(sample, "B_compactness")
        a, b = s["spearman_legacy"], s["spearman_amber"]
        if a != a or b != b:
            return float("nan")
        return abs(a) - abs(b)

    b2_point = _b2_stat(trials)
    b2_boot = _bootstrap_over_targets(trials, _b2_stat)
    b2_ci = [float(np.percentile(b2_boot, 2.5)), float(np.percentile(b2_boot, 97.5))] if len(
        b2_boot) else [float("nan")] * 2
    print(f"\nF-B2 COMPACTNESS SPECIALISATION  |sens_B(Legacy)| - |sens_B(AMBER)|"
          f" (Spearman rho vs dRg)")
    print(f"  point estimate = {b2_point:.4f}   bootstrap 95% CI = [{b2_ci[0]:.4f}, {b2_ci[1]:.4f}]")
    fb2_fires = bool(b2_ci[0] > 0.0)
    fb2_refuted = bool(b2_ci[1] < 0.0)
    print(f"  F-B2: {'SUPPORTED (Legacy > AMBER)' if fb2_fires else 'REFUTED (AMBER > Legacy)' if fb2_refuted else 'NOT MEASURED (CI spans 0)'}")

    # ---- per-target-permutation null, one representative axis pair, reported not gating ----
    rng = np.random.default_rng(0)

    def null_spearman_for(ax, energy_key):
        sub = [t for t in trials if t["axis"] == ax]
        proxy = PROXY_OF[ax]
        by_pdb = {}
        for t in sub:
            by_pdb.setdefault(t["pdb"], []).append(t)
        vals = []
        for _ in range(500):
            xs, ys = [], []
            for pdb, ts in by_pdb.items():
                x = [t[proxy] for t in ts]
                y = [t[energy_key] for t in ts]
                perm = rng.permutation(len(y))
                xs.extend(x); ys.extend([y[k] for k in perm])
            vals.append(_spearman(xs, ys))
        return np.asarray(vals, float)

    print("\nWITHIN-TARGET PERMUTATION NULL (proxy vs dE pairing shuffled within target, axis "
          "fixed), 500 draws, |rho| 95th percentile:")
    for ax in AXES:
        nl = null_spearman_for(ax, "zL")
        na = null_spearman_for(ax, "zA")
        print(f"  {ax:<20} null|rho| p95  Legacy {np.nanpercentile(np.abs(nl),95):.4f}   "
              f"AMBER {np.nanpercentile(np.abs(na),95):.4f}   (observed above)")

    result = {
        "source": os.path.basename(path), "n_targets": len(o["rows"]), "n_trials": len(trials),
        "sensitivity": sens,
        "F_B1": {"R_legacy": R_legacy, "R_amber": R_amber, "ratio": ratio, "ci95": ci,
                 "verdict": ("SUPPORTED" if fb1_fires else "REFUTED" if fb1_refuted else
                             "NOT_MEASURED")},
        "F_B2": {"point": b2_point, "ci95": b2_ci,
                 "verdict": ("SUPPORTED" if fb2_fires else "REFUTED" if fb2_refuted else
                             "NOT_MEASURED")},
    }
    outp = os.path.join(RESULTS, f"probe_perturb_analysis_{o['cfg_hash']}.json")
    with open(outp, "w") as fh:
        json.dump(result, fh, indent=1)
    print(f"\nwritten {outp}")
    return result


if __name__ == "__main__":
    a = sys.argv[1:]
    cmd = a[0] if a else "run"
    if cmd == "--smoke":
        run(subset=SUBSET[:2], out="_SMOKE_probe_perturb.json")
    elif cmd == "run":
        run()
    elif cmd == "analyse":
        analyse()
    else:
        raise SystemExit("usage: python -m s22.probe_perturb {run|analyse|--smoke}")
