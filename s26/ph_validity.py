"""s26/ph_validity.py -- THE HEAVY-ATOM VALIDITY AXIS OF THE PRODUCTION RELAXATION (lane PH).

Pre-registered in `s26/PREREG_validity_axis.md`. AMBER: one production relaxation per target,
gated bit-exact against the cached `amber_ca` / `amber_e1`. Native-free throughout (the panel
reads no native); folds are read for the CI only.

    python s26/ph_validity.py run [--n N]      # per-target cells, resumable
    python s26/ph_validity.py report
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

OUT = os.path.join(L.RESULTS, "ph_validity.json")
AXES = ("n_clash_2A", "n_clash_2p6A", "min_heavy", "bond_strain", "angle_strain", "rama_favoured",
        "rama_outlier", "cis_frac", "chirality_L_frac", "omega_dev")
KEYS = ("pdb", "n", "fold", "before", "after", "helix", "gate_dca", "gate_de", "gate_pass", "e0", "e1", "wall")
HELIX_PHI, HELIX_PSI = math.radians(-63.0), math.radians(-42.0)


def backbone_from_torsions(phi, psi):
    from core import geometry as geo
    bb = geo.build_backbone(np.asarray(phi, float), np.asarray(psi, float))
    return {k: np.asarray(v, float) for k, v in bb.items()}


def with_cb(bb: dict) -> dict:
    from core import geometry as geo
    out = {k: np.asarray(v, float) for k, v in bb.items() if k in ("N", "CA", "C", "O")}
    out["CB"] = np.array([geo.place_cb(out["N"][i], out["CA"][i], out["C"][i]) for i in range(len(out["CA"]))])
    return out


def panel_of(bb: dict, seq: str) -> dict:
    from s16 import energy_lib as EL
    p = EL.panel(bb, seq)
    return {k: float(p[k]) for k in AXES if k in p}


def run_target(t) -> dict:
    from core import amber as ar
    import torsion_lib2 as tl2
    rec = L.prod_record_nativefree(t["pdb"])
    seq = t["seq"]
    bb0 = backbone_from_torsions(rec["phi"], rec["psi"])
    before = panel_of(bb0, seq)
    helix = panel_of(backbone_from_torsions(np.full(t["n"], HELIX_PHI), np.full(t["n"], HELIX_PSI)), seq)
    tab = tl2.library_for(seq, 8, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    t0 = time.time()
    r = ar.refine_coords(seq, rep, bb0, k_restraint=10.0, steps=0, components=True, memo=False)
    wall = time.time() - t0
    ar.clear_cache()
    bb1 = with_cb({k: np.asarray(v, float) for k, v in r["backbone"].items()})
    after = panel_of(bb1, seq)
    dca = float(np.abs(np.asarray(r["ca"], float) - np.asarray(rec["amber_ca"], float)).max()) \
        if rec.get("amber_ca") is not None else float("nan")
    de = abs(float(r["energy"]) - float(rec["amber_e1"])) if rec.get("amber_e1") is not None else float("nan")
    return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "before": before, "after": after,
            "helix": helix, "gate_dca": dca, "gate_de": de, "gate_pass": bool(dca < 1e-6 and de < 1e-6),
            "e0": float(r["energy_initial"]), "e1": float(r["energy"]), "converged": bool(r["converged"]),
            "wall": round(wall, 2)}


def run(n: int = 0) -> None:
    tg = L.targets()
    if n:
        tg = tg[:n]
    done = L.read_cells("ph_validity")
    clk = L.Clock()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        r = run_target(t)
        L.write_cell("ph_validity", t["pdb"], r)
        print(f"[{k + 1}/{len(tg)}] {t['pdb']} gate dCA {r['gate_dca']:.1e} dE {r['gate_de']:.1e} "
              f"{'PASS' if r['gate_pass'] else 'FAIL'} clash2A {r['before']['n_clash_2A']:.0f}->{r['after']['n_clash_2A']:.0f} "
              f"bond {r['before']['bond_strain']:.4f}->{r['after']['bond_strain']:.4f} wall {r['wall']}s ({clk():.0f}s)", flush=True)
    rows = list(L.read_cells("ph_validity").values())
    L.save(OUT, {"what": "heavy-atom validity panel of the production emission before and after its relaxation, "
                         "with the constant alpha-helix reference; native-free", "rows": rows},
           rows=rows, complete_keys=KEYS, n_expected=126, module_file=__file__)


def report() -> dict:
    from s24 import stats_lib as ST
    rows = sorted(L.read_cells("ph_validity").values(), key=lambda r: r["pdb"])
    pdbs = [r["pdb"] for r in rows]; folds = L.folds_of(pdbs)
    out = {"n": len(rows), "gate_pass": int(sum(r["gate_pass"] for r in rows)),
           "gate_fail": [r["pdb"] for r in rows if not r["gate_pass"]],
           "gate_max_dca": float(max(r["gate_dca"] for r in rows)), "gate_max_de": float(max(r["gate_de"] for r in rows))}
    for ax in AXES:
        b = np.array([r["before"][ax] for r in rows]); a = np.array([r["after"][ax] for r in rows])
        h = np.array([r["helix"][ax] for r in rows])
        c = ST.compare(a, b, folds, names=pdbs, label=f"{ax}: relaxed minus built (negative = lower after)")
        print(ST.fmt(c))
        out[ax] = {"before": L.mean_se(b), "after": L.mean_se(a), "helix": L.mean_se(h), "compare": c}
    b2 = np.array([r["before"]["n_clash_2A"] for r in rows]); a2 = np.array([r["after"]["n_clash_2A"] for r in rows])
    out["targets_any_clash_2A_before"] = int((b2 > 0).sum()); out["targets_any_clash_2A_after"] = int((a2 > 0).sum())
    bs = np.array([r["after"]["bond_strain"] for r in rows])
    out["targets_bond_strain_after_over_0.05"] = int((bs > 0.05).sum())
    out["falsifier_clash_halved"] = bool(a2.mean() < 0.5 * b2.mean() and out["n_clash_2A"]["compare"]["ci95_fold"][1] < 0)
    print(json.dumps({k: v for k, v in out.items() if not isinstance(v, dict)}, indent=1))
    L.save(OUT, {"what": "heavy-atom validity panel, production emission before/after relaxation", "rows": rows,
                 "summary": out}, rows=rows, complete_keys=KEYS, n_expected=126, module_file=__file__)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("run", "report"))
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "run":
        run(a.n)
    else:
        report()
