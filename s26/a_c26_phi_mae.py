"""s26/a_c26_phi_mae.py -- lane A (Adversary): re-derive claim C26 from the S13 torsion cache.

C26 (`s26/EXAMINATION.md` section C; `docs/STATE_BRIEF_2026-09-12.md` section 5.6): phi MAE
36.1 deg for the full-sequence-context predictor (`p_grid`) against 36.4 deg for the
sequence-blind corpus marginal (`n_marg`); psi 62.4 against 72.8. The examination found no JSON
leaf for these; `s13/cache/tors_rows.npz` holds one JSON string per arm with 126 rows carrying
the per-residue absolute torsion errors `err_phi` and `err_psi`. This script pools them.

ORACLE DIAGNOSTIC: the cached errors were measured against the native torsions in S13. Nothing
here reads a structure, computes an RMSD, builds or selects anything. The record's statistic is
the MAE pooled over residues; the mean of per-target means is stored beside it.

    python s26/jobrun.py --agent A --tag CPU --name a_c26_phi_mae --est-ram 0.1 -- \
        python s26/a_c26_phi_mae.py

Writes `s26/results/a_c26_phi_mae.json` (`ST.save_atomic` provenance; complete gated on 126 rows).
"""
from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

CACHE = os.path.join(ROOT, "s13", "cache", "tors_rows.npz")
OUT = os.path.join(ROOT, "s26", "results", "a_c26_phi_mae.json")
KEY_ARMS = {"p_grid": "full sequence context + properties (the record's 36.1 / 62.4)",
            "n_marg": "sequence-blind corpus marginal (the record's 36.4 / 72.8)"}
RECORD = {"p_grid": (36.1, 62.4), "n_marg": (36.4, 72.8)}
NEED = ("pdb", "n", "fold", "n_phi", "n_psi", "p_grid_mae_phi", "n_marg_mae_phi",
        "p_grid_mae_psi", "n_marg_mae_psi")


def main():
    from s24 import stats_lib as ST
    z = np.load(CACHE, allow_pickle=True)
    arms = {a: json.loads(str(z[a])) for a in z.files}
    pdbs = [r["pdb"] for r in arms["p_grid"]]
    for a, rows in arms.items():
        assert [r["pdb"] for r in rows] == pdbs, a
    rows_out = []
    for i, pdb in enumerate(pdbs):
        base = arms["p_grid"][i]
        row = {"pdb": pdb, "n": int(base["n"]), "fold": int(base["fold"]),
               "n_phi": len(base["err_phi"]), "n_psi": len(base["err_psi"])}
        for a, rows in arms.items():
            r = rows[i]
            ep, es = np.abs(np.asarray(r["err_phi"], float)), np.abs(np.asarray(r["err_psi"], float))
            assert len(ep) == row["n_phi"] and len(es) == row["n_psi"], (a, pdb)
            row[f"{a}_sum_abs_phi"] = float(ep.sum()); row[f"{a}_sum_abs_psi"] = float(es.sum())
            row[f"{a}_mae_phi"] = float(ep.mean()); row[f"{a}_mae_psi"] = float(es.mean())
        rows_out.append(row)
    n_phi = sum(r["n_phi"] for r in rows_out); n_psi = sum(r["n_psi"] for r in rows_out)
    summary = {"n_targets": len(rows_out), "n_residues_phi": n_phi, "n_residues_psi": n_psi,
               "statistic": "MAE pooled over residues (sum of |err| / number of residues); "
                            "mean_of_target_means beside it", "arms": {}}
    for a in arms:
        pooled_phi = sum(r[f"{a}_sum_abs_phi"] for r in rows_out) / n_phi
        pooled_psi = sum(r[f"{a}_sum_abs_psi"] for r in rows_out) / n_psi
        summary["arms"][a] = {
            "pooled_mae_phi_deg": pooled_phi, "pooled_mae_psi_deg": pooled_psi,
            "mean_of_target_means_phi_deg": float(np.mean([r[f"{a}_mae_phi"] for r in rows_out])),
            "mean_of_target_means_psi_deg": float(np.mean([r[f"{a}_mae_psi"] for r in rows_out]))}
        if a in RECORD:
            rp, rs = RECORD[a]
            summary["arms"][a]["record_phi_psi"] = [rp, rs]
            summary["arms"][a]["matches_record_at_one_decimal"] = bool(
                round(pooled_phi, 1) == rp and round(pooled_psi, 1) == rs)
            summary["arms"][a]["role"] = KEY_ARMS[a]
    summary["claim_C26"] = {
        "phi_full_context_p_grid": summary["arms"]["p_grid"]["pooled_mae_phi_deg"],
        "phi_sequence_blind_n_marg": summary["arms"]["n_marg"]["pooled_mae_phi_deg"],
        "psi_full_context_p_grid": summary["arms"]["p_grid"]["pooled_mae_psi_deg"],
        "psi_sequence_blind_n_marg": summary["arms"]["n_marg"]["pooled_mae_psi_deg"],
        "record": "36.1 vs 36.4 deg phi, 62.4 vs 72.8 deg psi (s13/SPRINT13_DOSSIER.md:261-263)",
        "both_match_at_one_decimal": bool(summary["arms"]["p_grid"]["matches_record_at_one_decimal"]
                                          and summary["arms"]["n_marg"]["matches_record_at_one_decimal"])}
    for a, d in summary["arms"].items():
        print(f"{a:10} phi {d['pooled_mae_phi_deg']:7.3f}  psi {d['pooled_mae_psi_deg']:7.3f}  "
              f"(per-target means {d['mean_of_target_means_phi_deg']:7.3f} / {d['mean_of_target_means_psi_deg']:7.3f})"
              + (f"  record {d['record_phi_psi']} match={d['matches_record_at_one_decimal']}" if a in RECORD else ""))
    out = {"summary": summary, "rows": rows_out, "source": os.path.relpath(CACHE, ROOT),
           "label": "ORACLE DIAGNOSTIC (cached torsion errors against the native; no RMSD, no selection)"}
    ST.save_atomic(OUT, out, complete_keys=NEED, rows=rows_out, n_expected=126, module_file=__file__)
    print("wrote", os.path.relpath(OUT, ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
