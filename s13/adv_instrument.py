"""SPRINT 13 -- ADVERSARIAL AUDIT 1: the shared instrument.

Everything in this sprint stands on four objects:
    I.build_ca(phi, psi)                     -> the coordinates RMSD is measured on
    rep.build_coords(bitstring)              -> the coordinates core.amber / core.energy see
    tl2.library_for(seq, k, exclude_seq=seq)  -> the state table, claimed leakage-safe
    PerResidueTorsion._phi / _psi            -> the (n,k) tables every agent indexes directly

A defect in any one contaminates every claim at once.  Six independent checks, each of
which is a numerical measurement and not an argument.

    python -m s13.adv_instrument
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
import peptide_db as pdb                   # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402
import representations as reps             # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(RESULTS, exist_ok=True)
OUT = {}


def _wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


# --------------------------------------------------------------- A1 terminal torsions
def a1_terminal_torsions(tg, rng):
    """Are phi[0] and psi[n-1] free parameters that move nothing?"""
    rows = []
    for t in tg[:20]:
        n = t["n"]
        phi = rng.uniform(-np.pi, np.pi, n); psi = rng.uniform(-np.pi, np.pi, n)
        base = I.build_ca(phi, psi)
        d = {}
        for name, arr, idxs in (("phi0", "phi", [0]), ("psiLast", "psi", [n - 1]),
                                ("psi0", "psi", [0]), ("phiLast", "phi", [n - 1]),
                                ("phi_mid", "phi", [n // 2]), ("psi_mid", "psi", [n // 2])):
            p2, s2 = phi.copy(), psi.copy()
            (p2 if arr == "phi" else s2)[idxs] += 1.0     # 57 degrees
            d[name] = float(np.abs(I.build_ca(p2, s2) - base).max())
        rows.append(dict(pdb=t["pdb"], n=n, **d))
    agg = {k: float(np.max([r[k] for r in rows])) for k in rows[0] if k not in ("pdb", "n")}
    return {"per_target": rows, "max_abs_coord_change_over_20_targets": agg}


def a2_library_is_sequence_blind(tg, k=4):
    """How much sequence information does library_for(..., mode='class') actually carry?"""
    rows = []
    for t in tg:
        seq = t["seq"]
        tab = tl2.library_for(seq, k, seq)
        cls = reps.residue_classes(seq, len(seq))
        # distinct state-tables among residues
        flat = tab.reshape(len(seq), -1)
        uniq = np.unique(np.round(flat, 9), axis=0)
        rows.append({"pdb": t["pdb"], "n": t["n"], "n_distinct_state_tables": int(len(uniq)),
                     "n_distinct_classes": int(len(set(cls))),
                     "frac_GENERAL": float(np.mean([c == reps.CLASS_GENERAL for c in cls]))})
    return {"per_target": rows,
            "mean_distinct_state_tables": float(np.mean([r["n_distinct_state_tables"] for r in rows])),
            "mean_n_residues": float(np.mean([r["n"] for r in rows])),
            "equal_to_n_distinct_classes": int(sum(r["n_distinct_state_tables"] == r["n_distinct_classes"] for r in rows)),
            "mean_frac_GENERAL": float(np.mean([r["frac_GENERAL"] for r in rows]))}


def a3_holdout_bites(tg, k=4):
    """Does exclude_seq=seq move the table at all?  If not, 'leakage-safe' is decoration."""
    rows = []
    for t in tg[:40]:
        seq = t["seq"]
        held = tl2.library_for(seq, k, seq)
        full = tl2.library_for(seq, k, "")
        # states are unordered; match greedily by nearest
        dmax = 0.0
        for i in range(len(seq)):
            A = held[i]; B = full[i]
            D = np.abs(_wrap(A[:, None, 0] - B[None, :, 0])) + np.abs(_wrap(A[:, None, 1] - B[None, :, 1]))
            dmax = max(dmax, float(D.min(1).max()))
        n_db_full = len(pdb.load()); n_db_held = len(pdb.holdout(seq))
        rows.append({"pdb": t["pdb"], "max_matched_state_shift_deg": float(np.degrees(dmax)),
                     "db_entries_removed": int(n_db_full - n_db_held), "db_entries_full": int(n_db_full)})
    return {"per_target": rows,
            "median_shift_deg": float(np.median([r["max_matched_state_shift_deg"] for r in rows])),
            "max_shift_deg": float(np.max([r["max_matched_state_shift_deg"] for r in rows])),
            "median_entries_removed": float(np.median([r["db_entries_removed"] for r in rows]))}


def a4_state_indexing(tg, k=4, rng=None):
    """PerResidueTorsion._phi/_psi vs the table, and bitstring round-trip vs build_coords."""
    from core import geometry as geo
    rows = []
    for t in tg[:15]:
        seq, n = t["seq"], t["n"]
        tab = tl2.library_for(seq, k, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        e_phi = float(np.abs(np.asarray(rep._phi) - tab[:, :, 0]).max())
        e_psi = float(np.abs(np.asarray(rep._psi) - tab[:, :, 1]).max())
        s = rng.integers(0, k, n)
        bits = rep.bitstring_from_states(np.asarray(s, int))
        # 1) does the bitstring decode back to the same states?
        back = np.asarray(rep.state_indices(bits), int)
        # 2) do the coordinates the energy models see match I.build_ca on the same states?
        ca_rep = np.asarray(rep.build_coords(bits)["CA"], float)
        phi = tab[np.arange(n), s, 0]; psi = tab[np.arange(n), s, 1]
        ca_ins = I.build_ca(phi, psi)
        bb = geo.build_backbone(phi, psi)
        ca_bb = np.asarray(bb["CA"], float)
        row = {"pdb": t["pdb"], "err_phi_table": e_phi, "err_psi_table": e_psi,
               "bit_roundtrip_ok": (None if back is None else bool(np.array_equal(back, s))),
               "max_dev_build_ca_vs_build_backbone": float(np.abs(ca_ins - ca_bb).max()),
               "rmsd_build_ca_vs_build_backbone": I.ca_rmsd(ca_ins, ca_bb)}
        if ca_rep is not None:
            row["max_dev_build_ca_vs_rep_build_coords"] = float(np.abs(ca_ins - ca_rep).max())
            row["rmsd_build_ca_vs_rep_build_coords"] = I.ca_rmsd(ca_ins, ca_rep)
        rows.append(row)
    agg = {}
    for key in rows[0]:
        if key in ("pdb", "bit_roundtrip_ok"):
            continue
        agg["max_" + key] = float(np.max([r[key] for r in rows if key in r]))
    agg["bit_roundtrip_all_ok"] = all(r["bit_roundtrip_ok"] in (True, None) for r in rows)
    return {"per_target": rows, "agg": agg}


def a5_amber_sees_same_coords(tg, k=4, rng=None):
    """Does core.amber score the CA trace the RMSD is computed on?"""
    import core.amber as AM
    rows = []
    for t in tg[:4]:
        seq, n = t["seq"], t["n"]
        tab = tl2.library_for(seq, k, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        s = rng.integers(0, k, n)
        for _ in range(180):
            if AM.memory_percent() < AM.MEMORY_LIMIT_PERCENT - 2.0:
                break
            import time as _t; _t.sleep(10.0)
        r = AM.single_point(seq, rep, s)
        ca_amber = np.asarray(r["ca"], float)
        phi = tab[np.arange(n), s, 0]; psi = tab[np.arange(n), s, 1]
        ca_ins = I.build_ca(phi, psi)
        rows.append({"pdb": t["pdb"],
                     "rmsd_amber_ca_vs_instrument_ca": I.ca_rmsd(ca_amber, ca_ins),
                     "max_dev": float(np.abs(ca_amber - ca_ins).max())})
    return {"per_target": rows,
            "max_rmsd": float(np.max([r["rmsd_amber_ca_vs_instrument_ca"] for r in rows]))}


def a6_snap_uses_dead_angles(tg, k=4):
    """Does `ceiling.snap` choose terminal states on torsions the builder ignores?"""
    from s13 import ceiling as CE
    rows = []
    for t in tg[:40]:
        seq, n = t["seq"], t["n"]
        tab = tl2.library_for(seq, k, seq)
        PHI, PSI = tab[:, :, 0].copy(), tab[:, :, 1].copy()
        p = pdb.by_pdb(t["pdb"])
        phi0 = np.asarray(p.phi, float); psi0 = np.asarray(p.psi, float)
        s_full = CE.snap(PHI, PSI, phi0, psi0)
        # snap that ignores the dead angles at the termini
        d = np.abs(_wrap(PHI - phi0[:, None])) + np.abs(_wrap(PSI - psi0[:, None]))
        d0 = np.abs(_wrap(PSI[0] - psi0[0]))                        # residue 0: only psi lives
        dl = np.abs(_wrap(PHI[n - 1] - phi0[n - 1]))                # residue n-1: only phi lives
        d2 = d.copy(); d2[0] = d0; d2[n - 1] = dl
        s_live = np.argmin(d2, axis=1)
        u = I.load_univ(t["pdb"]); nat = u["nat_ca"]
        b = lambda s: I.build_ca(PHI[np.arange(n), s], PSI[np.arange(n), s])  # noqa: E731
        rows.append({"pdb": t["pdb"],
                     "snap_full": I.ca_rmsd(b(s_full), nat),
                     "snap_liveangles": I.ca_rmsd(b(s_live), nat),
                     "terminal_states_differ": int((s_full[0] != s_live[0]) + (s_full[n - 1] != s_live[n - 1]))})
    a = np.array([r["snap_full"] for r in rows]); c = np.array([r["snap_liveangles"] for r in rows])
    return {"per_target": rows, "mean_snap_full": float(a.mean()),
            "mean_snap_liveangles": float(c.mean()),
            "mean_diff": float((c - a).mean()),
            "n_targets_terminal_state_changed": int(sum(r["terminal_states_differ"] > 0 for r in rows))}


def main():
    rng = np.random.default_rng(20260905)
    tg = I.targets()
    print(f"{len(tg)} targets")
    OUT["A1_terminal_torsions"] = a1_terminal_torsions(tg, rng)
    print("A1", json.dumps(OUT["A1_terminal_torsions"]["max_abs_coord_change_over_20_targets"]))
    OUT["A2_library_sequence_blind"] = a2_library_is_sequence_blind(tg)
    print("A2", {k: v for k, v in OUT["A2_library_sequence_blind"].items() if k != "per_target"})
    OUT["A4_state_indexing"] = a4_state_indexing(tg, rng=rng)
    print("A4", OUT["A4_state_indexing"]["agg"])
    OUT["A6_snap_dead_angles"] = a6_snap_uses_dead_angles(tg)
    print("A6", {k: v for k, v in OUT["A6_snap_dead_angles"].items() if k != "per_target"})
    OUT["A3_holdout"] = a3_holdout_bites(tg)
    print("A3", {k: v for k, v in OUT["A3_holdout"].items() if k != "per_target"})
    try:
        OUT["A5_amber_coords"] = a5_amber_sees_same_coords(tg, rng=rng)
        print("A5", {k: v for k, v in OUT["A5_amber_coords"].items() if k != "per_target"})
    except Exception as e:                                              # noqa: BLE001
        OUT["A5_amber_coords"] = {"error": repr(e)}
        print("A5 FAILED", repr(e))
    with open(os.path.join(RESULTS, "adv_instrument.json"), "w") as fh:
        json.dump(OUT, fh, indent=1)
    print("wrote adv_instrument.json")


if __name__ == "__main__":
    main()
