"""s26/ph_relief.py -- ROTAMER RELIEF OF THE AMBER SINGLE POINT (lane PH, S26).

Pre-registered in `s26/PREREG_rotamer_relief.md`. Part B runs only if the tournament ranks
`s26/IDEA_rotamer_relief.md` as a survivor.

    python s26/ph_relief.py probe --pdb 1A13 --m 5     # AMBER: relieve m members of one target
    python s26/ph_relief.py run [--n N]                # AMBER: every target, top-75, resumable cells
    python s26/ph_relief.py report                     # gated: B(ii) ranking and B(iii) reject

The relieved energy: one greedy sweep over the residues that have a chi1 (all but G, A, P), the
option set {builder default, 60, 180, 300} degrees, the lowest single point kept per residue.
Single points are genuine ff14SB/GBn2 evaluations with NO minimisation, through the fast path
(`_heavy_positions` + `_assemble` + one `getState`) that `s20.qb2_lib.AmberSP` uses, gated
bit-exact against `core.amber.refine_coords(k_restraint=0, steps=-1)` and against the cache.
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
from s26 import ph_reject as PR                               # noqa: E402

OPTIONS = (60.0, 180.0, 300.0)
NO_CHI1 = ("G", "A", "P")
RUN_JSON = os.path.join(L.RESULTS, "ph_relief_run.json")
REPORT_JSON = os.path.join(L.RESULTS, "ph_relief_report.json")
RUN_KEYS = ("pdb", "n", "fold", "sub", "e_raw", "e_relief", "chi1", "n_sp", "wall", "gate_max_rel", "gate_n")


# ------------------------------------------------------------------ pure logic (tested)
def option_set(default: float) -> list:
    """The builder's default first, then the staggered values; duplicates (mod 360) removed."""
    out = [float(default)]
    for o in OPTIONS:
        if all(abs(((o - x) + 180.0) % 360.0 - 180.0) > 1e-6 for x in out):
            out.append(float(o))
    return out


def greedy_relief(energy_fn, seq: str, defaults: dict) -> tuple:
    """One greedy sweep. `energy_fn(chi1_dict) -> float`; `defaults[i]` = the builder's chi1 for
    residue i. Returns (E_relief, chi1 dict, n_evaluations). Never worse than the raw energy."""
    chi1 = {}
    e_cur = energy_fn(chi1)
    n_ev = 1
    for i, aa in enumerate(seq):
        if aa in NO_CHI1 or i not in defaults:
            continue
        best_e, best_o = e_cur, None
        for o in option_set(defaults[i]):
            trial = dict(chi1); trial[i] = o
            e = energy_fn(trial); n_ev += 1
            if e < best_e:
                best_e, best_o = e, o
        if best_o is not None:
            chi1[i] = best_o; e_cur = best_e
    return float(e_cur), chi1, n_ev


# ------------------------------------------------------------------ the single point
class SinglePoint:
    """Genuine ff14SB/GBn2 single point of a built backbone with an optional chi1 override,
    no minimisation. Same arithmetic as `s20.qb2_lib.AmberSP._e`, plus chi1."""

    def __init__(self, seq: str):
        from core import amber as am
        import torsion_lib2 as tl2
        self.am = am
        self.seq = seq
        #: k = 4, the representation `s24/cache_amber` was built with (`s24/d2_amberscore.py:114`,
        #: `s13.qarch_lib.Space(pdb, 4).rep`).  `builder_for` keys its OpenMM context on
        #: `rep.n_states` and calibrates hydrogen frames on the representation's reference
        #: states, so a different k gives a different (legitimate) single point; the probe's
        #: first gate fired on exactly that (rel 0.47 against the cache at k = 8, 0.0 against
        #: `refine_coords`).  PREREG_rotamer_relief addendum 1.
        tab = tl2.library_for(seq, 4, seq)
        self.rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        self.H = am.builder_for(seq, self.rep, "CPU", 1)
        self.n_calls = 0

    def defaults(self) -> dict:
        """The builder's chi1 per residue (degrees), for residues that have one."""
        out = {}
        for i, aa in enumerate(self.seq):
            if aa in NO_CHI1:
                continue
            key = self.am.THREE.get(aa, "GLY")
            chis = self.am.CHI_ANGLES.get(key, ())
            if chis:
                out[i] = float(chis[0])
        return out

    def energy(self, coords: dict, chi1: dict | None) -> float:
        from openmm import unit
        H = self.H
        pos = H._assemble(H._heavy_positions(coords, chi1=chi1 or None))
        ctx = H.context
        ctx.setPositions(pos * unit.nanometer)
        ctx.setParameter("k_rest", 0.0)
        self.n_calls += 1
        return float(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))

    def verify(self, coords: dict, chi1: dict | None = None) -> float:
        """Relative difference against `refine_coords(k=0, steps=-1)`; 0.0 is required."""
        ref = float(self.am.refine_coords(self.seq, self.rep, coords, k_restraint=0.0, steps=-1,
                                          tolerance=1e9, chi1=chi1, threads=1, memo=False)["energy"])
        got = self.energy(coords, chi1)
        return abs(ref - got) / max(1.0, abs(ref))


# ------------------------------------------------------------------ one target
def relieve_target(t, m: int = 0) -> dict:
    from core import geometry as geo
    pool = PR.load_pool(t["pdb"], oracle=False)
    u = L.univ_nativefree(t["pdb"], keys=("PHI", "PSI", "order"))
    p = L.pool_idx_from_order(u)
    PHI, PSI = u["PHI"][p], u["PSI"][p]
    sub = pool["sub"] if not m else pool["sub"][:m]
    sp = SinglePoint(pool["seq"])
    defaults = sp.defaults()
    e_raw, e_rel, chis, nsp, walls = [], [], [], [], []
    gate = []
    for k, j in enumerate(sub):
        bb = geo.build_backbone(PHI[j], PSI[j])
        if k < 4:
            # G1: the fast path must reproduce refine_coords AND the cache bit-for-bit
            rel = sp.verify(bb, None)
            cached = float(pool["e"][j])
            rel2 = abs(sp.energy(bb, None) - cached) / max(1.0, abs(cached))
            gate.append(max(rel, rel2))
            assert max(rel, rel2) == 0.0, f"{t['pdb']} member {j}: single point not bit-exact ({rel:.2e}, {rel2:.2e})"
        t0 = time.time()
        e0 = float(pool["e"][j])
        er, chi1, n_ev = greedy_relief(lambda c: sp.energy(bb, c), pool["seq"], defaults)
        assert er <= e0 + 1e-9 * max(1.0, abs(e0)), f"{t['pdb']} member {j}: relief above raw"
        e_raw.append(e0); e_rel.append(er); chis.append({str(i): v for i, v in chi1.items()})
        nsp.append(n_ev); walls.append(round(time.time() - t0, 3))
    return {"pdb": t["pdb"], "n": int(t["n"]), "fold": int(t["fold"]), "sub": [int(x) for x in sub],
            "e_raw": e_raw, "e_relief": e_rel, "chi1": chis, "n_sp": nsp, "wall": walls,
            "gate_max_rel": float(max(gate)), "gate_n": len(gate), "n_chi1_residues": len(defaults),
            "frac_over_1e4_raw": float(np.mean(np.array(e_raw) > 1e4)),
            "frac_over_1e4_relief": float(np.mean(np.array(e_rel) > 1e4)),
            "frac_over_1e6_raw": float(np.mean(np.array(e_raw) > 1e6)),
            "frac_over_1e6_relief": float(np.mean(np.array(e_rel) > 1e6))}


def probe(pdb: str, m: int = 5) -> dict:
    t = {x["pdb"]: x for x in L.targets()}[pdb]
    r = relieve_target(t, m=m)
    print(json.dumps({"pdb": pdb, "m": m, "gate_max_rel": r["gate_max_rel"], "gate_n": r["gate_n"],
                      "n_chi1_residues": r["n_chi1_residues"], "n_sp_per_member": r["n_sp"],
                      "wall_per_member_s": r["wall"], "e_raw": [round(x, 1) for x in r["e_raw"]],
                      "e_relief": [round(x, 1) for x in r["e_relief"]],
                      "extrapolated_126x75_hours": round(126 * 75 * float(np.mean(r["wall"])) / 3600.0, 2)}, indent=1))
    return r


def run(n: int = 0) -> None:
    tg = L.targets()
    if n:
        tg = tg[:n]
    done = L.read_cells("ph_relief")
    clk = L.Clock()
    for k, t in enumerate(tg):
        if t["pdb"] in done:
            continue
        r = relieve_target(t)
        L.write_cell("ph_relief", t["pdb"], r)
        print(f"[{k + 1}/{len(tg)}] {t['pdb']} over1e4 raw {r['frac_over_1e4_raw']:.2f} -> relief "
              f"{r['frac_over_1e4_relief']:.2f}  sp/member {np.mean(r['n_sp']):.0f}  wall {sum(r['wall']):.0f}s ({clk():.0f}s)",
              flush=True)
    rows = list(L.read_cells("ph_relief").values())
    summ = {"n": len(rows),
            "frac_over_1e4_raw": L.mean_se([r["frac_over_1e4_raw"] for r in rows]),
            "frac_over_1e4_relief": L.mean_se([r["frac_over_1e4_relief"] for r in rows]),
            "frac_over_1e6_raw": L.mean_se([r["frac_over_1e6_raw"] for r in rows]),
            "frac_over_1e6_relief": L.mean_se([r["frac_over_1e6_relief"] for r in rows]),
            "gate_max_rel": float(max(r["gate_max_rel"] for r in rows)),
            "falsifier_Bi_fires": bool(np.mean([r["frac_over_1e4_relief"] for r in rows])
                                       >= 0.5 * np.mean([r["frac_over_1e4_raw"] for r in rows]))}
    L.save(RUN_JSON, {"what": "rotamer relief of the AMBER single point on the shipped top-75 (native-free)",
                      "rows": rows, "summary": summ}, rows=rows, complete_keys=RUN_KEYS, n_expected=126,
           module_file=__file__)
    print(json.dumps(summ, indent=1))


# ------------------------------------------------------------------ report (gated)
def report() -> dict:
    L.require_gate("ph_relief report")
    from s24 import stats_lib as ST
    from scipy.stats import spearmanr
    rows = sorted(L.read_cells("ph_relief").values(), key=lambda r: r["pdb"])
    pdbs = [r["pdb"] for r in rows]; folds = L.folds_of(pdbs)
    rho_raw, rho_rel, rho_raw_ib, rho_rel_ib = [], [], [], []
    for r in rows:
        u = L.univ_oracle(r["pdb"]); p = L.pool_idx_from_order(u)
        sub = np.asarray(r["sub"], int)
        rr = u["rr"][p][sub]
        er, el = np.array(r["e_raw"]), np.array(r["e_relief"])
        rho_raw.append(spearmanr(er, rr)[0]); rho_rel.append(spearmanr(el, rr)[0])
        ib = rr <= rr.min() + 3.0
        if ib.sum() >= 5:
            rho_raw_ib.append(spearmanr(er[ib], rr[ib])[0]); rho_rel_ib.append(spearmanr(el[ib], rr[ib])[0])
        else:
            rho_raw_ib.append(np.nan); rho_rel_ib.append(np.nan)
    out = {"rho_raw": L.mean_se(rho_raw), "rho_relief": L.mean_se(rho_rel),
           "rho_raw_inband": L.mean_se(rho_raw_ib), "rho_relief_inband": L.mean_se(rho_rel_ib)}
    ok = np.isfinite(rho_rel_ib) & np.isfinite(rho_raw_ib)
    for lab, a, b, m in (("rho(E_relief, ORACLE d) minus rho(E_raw, ORACLE d), whole top-75", np.array(rho_rel), np.array(rho_raw), np.ones(len(rows), bool)),
                         ("in-band rho(E_relief) minus rho(E_raw)", np.array(rho_rel_ib), np.array(rho_raw_ib), ok)):
        c = ST.compare(a[m], b[m], folds[m], names=[p for p, o in zip(pdbs, m) if o], label=lab)
        print(ST.fmt(c)); out[lab] = c
    L.save(REPORT_JSON, {"what": "B(ii): does the relieved single point rank? (ORACLE evaluation)", "report": out},
           module_file=__file__)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("probe", "run", "report"))
    ap.add_argument("--pdb", default="1A13")
    ap.add_argument("--m", type=int, default=5)
    ap.add_argument("--n", type=int, default=0)
    a = ap.parse_args()
    if a.mode == "probe":
        probe(a.pdb, a.m)
    elif a.mode == "run":
        run(a.n)
    else:
        report()
