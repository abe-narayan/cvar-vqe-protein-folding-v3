"""s26/ph_c3.py -- C3, THE MATCHED-RANDOM CONTROL FOR AMBER REFINEMENT (lane PH, S26).

Pre-registered in `s26/PREREG_c3_control.md`.  Proposal C, Part 3.3 of the campaign prompt.

WHAT IS MEASURED.  The production relaxation (`core.pipeline.relax`: `refine_coords(k_restraint=
10, steps=0)`, converged) moves the built chain `ca` (basis `rmsd_arm`, 3.2148) to `amber_ca`
(basis `rmsd_full`, 3.2355), +0.0207 A [+0.0143, +0.0276] on the 126 dev targets (state brief
section 4).  S16 measured, on a DIFFERENT input (the all-atom coordinate average, k=30), that a
random displacement of matched magnitude is at least as accurate as AMBER's and that AMBER's
displacement points slightly away from the truth (cos -0.052).  C3 asks the same question on the
PRODUCTION input and the PRODUCTION operator, from the persisted coordinates, with no new AMBER
compute (stage 1), and later on lane P's best C2 rung output (stage 2, AMBER compute).

THE CONTROLS, matched in the operator's space (a displacement's control moves the same distance):
    rand      S16's construction reproduced (`s16/repair_report.py:228-238`): an isotropic Gaussian
              direction on the CA trace, the six rigid-body components removed by projection onto
              `s15.align_lib.rigid_basis`, scaled to the SAME per-atom RMS displacement AMBER
              applied (measured after superposing `amber_ca` onto `ca`).  16 draws, mean.
    member    S20's `c_land_null` form in coordinate space: the same RMS magnitude along the
              straight line from `ca` toward another randomly chosen member of the shipped K=500
              pool (superposed onto `ca`).  Realisable, native-free, zero information.  16 draws.
    zero-information anchor: `ca` itself (do nothing).
DIFFERENCES FROM S16, stated: (i) input is the built chain, not the all-atom average; (ii) the
displacement magnitude is measured on CA after superposition, which is the quantity the control
reproduces (S16 used `moved_ca`, the same definition); (iii) 16 draws instead of 3.

ENDPOINTS (gated): full-chain CA-RMSD to the native of `ca + g` for every draw; AMBER minus each
control paired per target, `ST.compare`, fold CI beside iid, MDE, W/L, concentration null; the
ORACLE cosine of AMBER's displacement with the true residual beside the controls' cosines.

THE VALIDITY AXIS (native-free): from the production cache, the built chain's own AMBER energy
`amber_e0` (its strain before relaxation), the relaxed energy `amber_e1`, `amber_strain_after`
(bond+angle after), `amber_moved` (restraint RMSD, N/CA/C), and a CA-level axis (virtual bond
statistics and CA contacts) for `ca` and `amber_ca`.  The heavy-atom axis (clash count, bond
deviation) needs the relaxed backbone, which the cache does not hold; stage 2 measures it.

DECISION RULE (from the prompt): if the relaxation's gain over the random control is below its MDE,
the presentation drops "refine with physics" as an accuracy step and keeps it as a validity step.

    python s26/ph_c3.py nativefree                 # before the gate: magnitudes + validity scalars
    python s26/ph_c3.py stage1                     # gated: the controls on the production chain
    python s26/ph_c3.py probe --input F --pdb P    # stage 2 probe, one target, under jobrun AMBER
    python s26/ph_c3.py stage2 --input F           # stage 2, resumable per target, AMBER tag
"""
from __future__ import annotations

import argparse
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

from s26 import ph_lib as L                                   # noqa: E402

N_DRAWS = 16
NF_JSON = os.path.join(L.RESULTS, "ph_c3_nativefree.json")
S1_JSON = os.path.join(L.RESULTS, "ph_c3_stage1.json")
S2_JSON = os.path.join(L.RESULTS, "ph_c3_stage2.json")
AMBER_K = 10.0          # core.pipeline.Config.amber_k
AMBER_STEPS = 0         # core.pipeline.Config.amber_steps (0 = converged)
TORSION_WINDOW = 8      # core.pipeline.Config.torsion_window


# ------------------------------------------------------------------ native-free part
NF_KEYS = ("pdb", "n", "fold", "mag_sup", "mag_raw", "amber_moved", "e0", "e1", "converged",
           "strain_after", "ca_bond_mean", "ca_bond_min", "ca_bond_max", "amb_bond_mean",
           "amb_bond_min", "amb_bond_max", "ca_contacts_4A", "amb_contacts_4A", "rg_ca", "rg_amb")


def displacement(ca: np.ndarray, out: np.ndarray):
    """v = (operator output superposed onto its input) - input; returns (v, per-atom RMS)."""
    ca = np.asarray(ca, float); out = np.asarray(out, float)
    v = L.superpose_onto(out, ca) - ca
    return v, float(np.sqrt((v ** 2).sum(1).mean()))


def nativefree_target(t) -> dict:
    rec = L.prod_record_nativefree(t["pdb"])
    ca = np.asarray(rec["ca"], float)
    if rec.get("amber_ca") is None:
        return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "amber_err": rec.get("amber_err"),
                "mag_sup": float("nan")}
    amb = np.asarray(rec["amber_ca"], float)
    v, mag = displacement(ca, amb)
    raw = float(np.sqrt(((amb - ca) ** 2).sum(1).mean()))
    b0, b1 = L.consecutive_ca(ca), L.consecutive_ca(amb)
    e1 = float(rec["amber_e1"])
    return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]),
            "mag_sup": mag, "mag_raw": raw, "amber_moved": float(rec["amber_moved"]),
            "e0": float(rec["amber_e0"]), "e1": e1,
            "converged": bool(np.isfinite(e1) and e1 <= 1000.0),
            "strain_after": float(rec["amber_strain_after"]),
            "ca_bond_mean": float(b0.mean()), "ca_bond_min": float(b0.min()), "ca_bond_max": float(b0.max()),
            "amb_bond_mean": float(b1.mean()), "amb_bond_min": float(b1.min()), "amb_bond_max": float(b1.max()),
            "ca_contacts_4A": int(L.ca_contacts(ca, 4.0)), "amb_contacts_4A": int(L.ca_contacts(amb, 4.0)),
            "ca_contacts_3p8A": int(L.ca_contacts(ca, 3.8)), "amb_contacts_3p8A": int(L.ca_contacts(amb, 3.8)),
            "rg_ca": float(L.rg_of(ca)), "rg_amb": float(L.rg_of(amb)),
            "amber_err": rec.get("amber_err")}


def nativefree() -> dict:
    tg = L.targets()
    rows = [nativefree_target(t) for t in tg]
    ok = [r for r in rows if np.isfinite(r.get("mag_sup", float("nan")))]
    summ = {"n": len(rows), "n_with_amber": len(ok),
            "mag_sup": L.mean_se([r["mag_sup"] for r in ok]),
            "mag_raw": L.mean_se([r["mag_raw"] for r in ok]),
            "amber_moved": L.mean_se([r["amber_moved"] for r in ok]),
            "e0": L.mean_se([r["e0"] for r in ok]), "e0_median": float(np.median([r["e0"] for r in ok])),
            "e0_min": float(np.min([r["e0"] for r in ok])), "e0_max": float(np.max([r["e0"] for r in ok])),
            "e1": L.mean_se([r["e1"] for r in ok]), "e1_max": float(np.max([r["e1"] for r in ok])),
            "n_converged_le_1000": int(sum(r["converged"] for r in ok)),
            "strain_after": L.mean_se([r["strain_after"] for r in ok]),
            "ca_bond_mean": L.mean_se([r["ca_bond_mean"] for r in ok]),
            "amb_bond_mean": L.mean_se([r["amb_bond_mean"] for r in ok]),
            "amb_bond_min": float(np.min([r["amb_bond_min"] for r in ok])),
            "amb_bond_max": float(np.max([r["amb_bond_max"] for r in ok])),
            "ca_contacts_4A": L.mean_se([r["ca_contacts_4A"] for r in ok]),
            "amb_contacts_4A": L.mean_se([r["amb_contacts_4A"] for r in ok]),
            "d_rg_amb_minus_ca": L.mean_se([r["rg_amb"] - r["rg_ca"] for r in ok]),
            "frac_e0_over_1e4": float(np.mean([r["e0"] > 1e4 for r in ok]))}
    L.save(NF_JSON, {"what": "C3 native-free part: production AMBER displacement magnitudes and "
                             "the validity scalars from bench_results/cache/1fc9f2dcf489e2fb",
                     "rows": rows, "summary": summ},
           rows=rows, complete_keys=NF_KEYS, n_expected=126, module_file=__file__)
    print(json.dumps(summ, indent=1))
    return summ


# ------------------------------------------------------------------ stage 1 (gated)
S1_KEYS = ("pdb", "n", "fold", "rmsd_arm", "rmsd_full", "mag", "cos_amber", "rand_rmsd",
           "rand_cos", "member_rmsd", "member_cos", "pred_orthogonal")


def controls_for(ca: np.ndarray, out: np.ndarray, nat: np.ndarray, Wpool: np.ndarray,
                 pdb: str, n_draws: int = N_DRAWS) -> dict:
    """ORACLE DIAGNOSTIC (reads `nat`). The operator's displacement, its cosine with the true
    residual, and the two matched-magnitude controls, all on one input `ca`."""
    ca = np.asarray(ca, float); n = len(ca)
    v, mag = displacement(ca, out)
    r = L.superpose_onto(nat, ca) - ca
    rn = np.linalg.norm(r)

    def cos(a):
        na = np.linalg.norm(a)
        return float((a * r).sum() / (na * rn)) if na > 0 and rn > 0 else float("nan")

    d_out = L.ca_rmsd(out, nat)
    rng = L.stable_rng(pdb, "c3rand")
    rand_r, rand_c = [], []
    for _ in range(n_draws):
        g = L.random_displacement(ca, mag, rng)
        rand_r.append(L.ca_rmsd(ca + g, nat)); rand_c.append(cos(g))
    rng2 = L.stable_rng(pdb, "c3member")
    mem_r, mem_c, over = [], [], 0
    for _ in range(n_draws):
        j = int(rng2.integers(0, len(Wpool)))
        w = L.superpose_onto(Wpool[j], ca)
        d = w - ca
        dn = np.linalg.norm(d)
        s = mag * math.sqrt(n) / max(dn, 1e-12)
        over += int(s > 1.0)
        g = d * s
        mem_r.append(L.ca_rmsd(ca + g, nat)); mem_c.append(cos(g))
    # the exact identity n*RMSD_after^2 = |r|^2 - 2 v.r + |v|^2 (S16 section 3.4)
    pred = math.sqrt(max((rn ** 2 - 2 * (v * r).sum() + (v ** 2).sum()) / n, 0.0))
    pred_orth = math.sqrt(max((rn ** 2 + (v ** 2).sum()) / n, 0.0))
    return {"mag": mag, "cos_amber": cos(v), "rmsd_out": d_out, "pred_identity": pred,
            "pred_orthogonal": pred_orth,
            "rand_rmsd": float(np.mean(rand_r)), "rand_rmsd_sd": float(np.std(rand_r, ddof=1)),
            "rand_cos": float(np.mean(rand_c)),
            "member_rmsd": float(np.mean(mem_r)), "member_rmsd_sd": float(np.std(mem_r, ddof=1)),
            "member_cos": float(np.mean(mem_c)), "member_overshoot_draws": over,
            "n_draws": n_draws}


def stage1_target(t) -> dict:
    L.require_gate("ph_c3 stage1")
    from s12 import instrument as I
    rec = L.prod_record_oracle(t["pdb"])
    u = L.univ_oracle(t["pdb"])
    p = L.pool_idx_from_order(u)
    nat = u["nat_ca"]
    ca = np.asarray(rec["ca"], float)
    if rec.get("amber_ca") is None:
        return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "skipped": rec.get("amber_err")}
    amb = np.asarray(rec["amber_ca"], float)
    ra, rf = I.ca_rmsd(ca, nat), I.ca_rmsd(amb, nat)
    assert abs(ra - rec["rmsd_arm"]) < 1e-6 and abs(rf - rec["rmsd_full"]) < 1e-6, t["pdb"]
    c = controls_for(ca, amb, nat, u["W"][p], t["pdb"])
    return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "basis_in": "built_chain",
            "basis_out": "relaxed_chain", "rmsd_arm": float(ra), "rmsd_full": float(rf), **c}


def _stats(rows, label_prefix, out):
    from s24 import stats_lib as ST
    rows = [r for r in rows if "skipped" not in r]
    pdbs = [r["pdb"] for r in rows]; folds = L.folds_of(pdbs)
    A = np.array([r["rmsd_full"] for r in rows]); B = np.array([r["rmsd_arm"] for r in rows])
    R = np.array([r["rand_rmsd"] for r in rows]); Mm = np.array([r["member_rmsd"] for r in rows])
    cA = np.array([r["cos_amber"] for r in rows]); cR = np.array([r["rand_cos"] for r in rows])
    cM = np.array([r["member_cos"] for r in rows])
    for lab, a, b in ((f"{label_prefix} AMBER (relaxed) minus do-nothing (built chain)", A, B),
                      (f"{label_prefix} AMBER minus matched-magnitude RANDOM (16 draws)", A, R),
                      (f"{label_prefix} AMBER minus matched-magnitude TOWARD-MEMBER (16 draws)", A, Mm),
                      (f"{label_prefix} RANDOM minus do-nothing", R, B),
                      (f"{label_prefix} TOWARD-MEMBER minus do-nothing", Mm, B),
                      (f"{label_prefix} ORACLE cos: AMBER minus RANDOM", cA, cR),
                      (f"{label_prefix} ORACLE cos: AMBER minus TOWARD-MEMBER", cA, cM)):
        c = ST.compare(a, b, folds, names=pdbs, label=lab)
        print(ST.fmt(c))
        out[lab] = c
    out["cos_amber"] = L.mean_se(cA); out["cos_rand"] = L.mean_se(cR); out["cos_member"] = L.mean_se(cM)
    out["frac_cos_amber_positive"] = float((cA > 0).mean())
    out["mag"] = L.mean_se([r["mag"] for r in rows])
    out["identity_check_max_abs"] = float(np.max([abs(r["pred_identity"] - r["rmsd_full"]) for r in rows]))
    out["orthogonal_move_cost"] = L.mean_se([r["pred_orthogonal"] - r["rmsd_arm"] for r in rows])
    return out


def stage1() -> dict:
    L.require_gate("ph_c3 stage1")
    tg = L.targets()
    rows, clk = [], L.Clock()
    for k, t in enumerate(tg):
        r = stage1_target(t)
        rows.append(r)
        if "skipped" not in r:
            print(f"[{k + 1}/{len(tg)}] {r['pdb']} arm {r['rmsd_arm']:.3f} full {r['rmsd_full']:.3f} "
                  f"mag {r['mag']:.3f} cos {r['cos_amber']:+.3f} rand {r['rand_rmsd']:.3f} "
                  f"member {r['member_rmsd']:.3f} ({clk():.0f}s)", flush=True)
    summ = _stats(rows, "stage1", {})
    L.save(S1_JSON, {"what": "C3 stage 1: AMBER refinement vs matched-magnitude controls on the "
                             "production built chain; ORACLE evaluation of native-free operators",
                     "rows": rows, "summary": summ},
           rows=rows, complete_keys=S1_KEYS, n_expected=126, module_file=__file__)
    return summ


# ------------------------------------------------------------------ stage 2 (AMBER)
def relax_chain(seq: str, phi: np.ndarray, psi: np.ndarray) -> dict:
    """The production relaxation, `core.pipeline._relax_inner` semantics: ideal-geometry
    backbone from (phi, psi), the target's own sequence held out of the torsion library,
    `refine_coords(k_restraint=10, steps=0, components=True)`.  Returns CA, backbone, energies,
    the convergence flag and `s16.energy_lib.panel` on the relaxed backbone."""
    from core import geometry as geo
    from core import amber as ar
    import torsion_lib2 as tl2
    from s16 import energy_lib as EL
    BB = geo.build_backbone_batch(np.asarray(phi, float)[None], np.asarray(psi, float)[None])
    cd = {k: v[0] for k, v in BB.items()}
    tab = tl2.library_for(seq, TORSION_WINDOW, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    t0 = time.time()
    r = ar.refine_coords(seq, rep, cd, k_restraint=AMBER_K, steps=AMBER_STEPS, components=True,
                         memo=False)
    wall = time.time() - t0
    bb = {k: np.asarray(v, float) for k, v in r["backbone"].items()}
    pan_in = EL.panel(cd, seq)
    pan_out = EL.panel(bb, seq)
    ar.clear_cache()
    return {"ca": np.asarray(r["ca"], float), "backbone": bb, "e0": float(r["energy_initial"]),
            "e1": float(r["energy"]), "moved": float(r["restraint_rmsd"]),
            "converged": bool(r["converged"]), "converge_reason": r.get("converge_reason", ""),
            "components": r.get("components"), "wall": wall,
            "panel_in": pan_in, "panel_out": pan_out}


def _load_input(path: str) -> dict:
    """Lane P's output: a JSON with rows carrying `pdb`, `phi`, `psi` (radians, length n)."""
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    rows = d["rows"] if isinstance(d, dict) else d
    return {r["pdb"]: r for r in rows}


def probe(path: str, pdb: str) -> dict:
    L.require_gate("ph_c3 probe")
    inp = _load_input(path)
    t = {x["pdb"]: x for x in L.targets()}[pdb]
    r = inp[pdb]
    out = relax_chain(t["seq"], np.asarray(r["phi"], float), np.asarray(r["psi"], float))
    print(json.dumps({"pdb": pdb, "e0": out["e0"], "e1": out["e1"], "moved": out["moved"],
                      "converged": out["converged"], "wall_s": round(out["wall"], 2),
                      "clash_in": out["panel_in"]["n_clash_2A"], "clash_out": out["panel_out"]["n_clash_2A"],
                      "bond_strain_out": out["panel_out"]["bond_strain"]}, indent=1))
    return out


def stage2(path: str, n: int = 0) -> None:
    L.require_gate("ph_c3 stage2")
    from s12 import instrument as I
    inp = _load_input(path)
    tg = [t for t in L.targets() if t["pdb"] in inp]
    if n:
        tg = tg[:n]
    done = L.read_cells("ph_c3")
    clk = L.Clock()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        r = inp[t["pdb"]]
        u = L.univ_oracle(t["pdb"]); p = L.pool_idx_from_order(u); nat = u["nat_ca"]
        out = relax_chain(t["seq"], np.asarray(r["phi"], float), np.asarray(r["psi"], float))
        ca_in = np.asarray(I.build_ca(np.asarray(r["phi"], float), np.asarray(r["psi"], float)), float)
        c = controls_for(ca_in, out["ca"], nat, u["W"][p], t["pdb"])
        row = {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "basis_in": "built_chain(P)",
               "basis_out": "relaxed_chain", "rmsd_arm": L.ca_rmsd(ca_in, nat),
               "rmsd_full": L.ca_rmsd(out["ca"], nat), "e0": out["e0"], "e1": out["e1"],
               "moved": out["moved"], "converged": out["converged"], "wall": out["wall"],
               "panel_in": out["panel_in"], "panel_out": out["panel_out"],
               "ca_out": out["ca"].tolist(), **c}
        L.write_cell("ph_c3", t["pdb"], row)
        print(f"[{k + 1}/{len(tg)}] {t['pdb']} arm {row['rmsd_arm']:.3f} full {row['rmsd_full']:.3f} "
              f"mag {row['mag']:.3f} rand {row['rand_rmsd']:.3f} member {row['member_rmsd']:.3f} "
              f"e {row['e0']:.3g}->{row['e1']:.1f} wall {row['wall']:.1f}s ({clk():.0f}s)", flush=True)
    rows = list(L.read_cells("ph_c3").values())
    summ = _stats(rows, "stage2", {}) if len(rows) >= 10 else {}
    L.save(S2_JSON, {"what": "C3 stage 2: production relaxation on lane P's best C2 rung output "
                             "with the matched-magnitude controls and the heavy-atom validity axis",
                     "input": path, "rows": rows, "summary": summ},
           rows=rows, complete_keys=S1_KEYS, n_expected=len(inp), module_file=__file__)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("nativefree", "stage1", "probe", "stage2"))
    ap.add_argument("--input", default=None)
    ap.add_argument("--pdb", default=None)
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "nativefree":
        nativefree()
    elif a.mode == "stage1":
        stage1()
    elif a.mode == "probe":
        probe(a.input, a.pdb)
    else:
        stage2(a.input, a.n)
