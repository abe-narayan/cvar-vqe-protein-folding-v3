"""s26/ph_relief_reject.py -- rotamer_relief B(iii): arm S (reject from the shipped top-75, no refill)
on the RELIEVED single point, point-cloud basis, with the matched-count random and permuted
controls of `s26/ph_reject.py`. Pre-registered in `s26/PREREG_rotamer_relief.md` section 5.
Gated (reads nat_ca to score). Arm R needs relieved energies beyond the top-75 and is not run
(section 5: not without a top-150 budget).

    python s26/ph_relief_reject.py
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s26 import ph_lib as L                                   # noqa: E402
from s26 import ph_reject as PR                               # noqa: E402

OUT = os.path.join(L.RESULTS, "ph_relief_reject.json")
N_DRAWS = 16


def main() -> dict:
    L.require_gate("ph_relief_reject")
    from s24 import stats_lib as ST
    cells = L.read_cells("ph_relief")
    rows = []
    for t in L.targets():
        pdb = t["pdb"]
        pool = PR.load_pool(pdb, oracle=True)
        W, nat, sub = pool["W"], pool["nat_ca"], pool["sub"]
        c = cells[pdb]
        assert [int(x) for x in c["sub"]] == [int(x) for x in sub], pdb
        e_rel = np.asarray(c["e_relief"], float)      # aligned to sub
        e_raw = np.asarray(c["e_raw"], float)
        anchor = PR.cloud_rmsd(W, sub, nat)
        row = {"pdb": pdb, "n": int(t["n"]), "fold": int(t["fold"]), "anchor": anchor, "arms": {}}
        for T in PR.THRESHOLDS:
            k = PR.tkey(T)
            keep_rel = sub[e_rel <= T]; keep_raw = sub[e_raw <= T]
            r_rel = int(len(sub) - len(keep_rel)); r_raw = int(len(sub) - len(keep_raw))
            S_rel = keep_rel if len(keep_rel) else sub.copy()
            S_raw = keep_raw if len(keep_raw) else sub.copy()
            rs = L.stable_rng(pdb, "relief_rands", k); rp = L.stable_rng(pdb, "relief_perm", k)
            rand = [PR.random_shrink(sub, r_rel, rs) for _ in range(N_DRAWS)]
            rand = [x if len(x) else sub.copy() for x in rand]
            perm = []
            for _ in range(N_DRAWS):
                ep = e_rel[rp.permutation(len(e_rel))]
                kp = sub[ep <= T]; perm.append(kp if len(kp) else sub.copy())
            row["arms"][k] = {"n_reject_relief": r_rel, "n_reject_raw": r_raw,
                              "S_relief": PR.cloud_rmsd(W, S_rel, nat), "S_raw": PR.cloud_rmsd(W, S_raw, nat),
                              "RANDS": float(np.mean([PR.cloud_rmsd(W, x, nat) for x in rand])),
                              "PERMS": float(np.mean([PR.cloud_rmsd(W, x, nat) for x in perm])),
                              "S_empty": bool(len(keep_rel) == 0)}
        rows.append(row)
        print(f"{pdb} anchor {anchor:.3f} " + " ".join(f"{k}:S_rel {row['arms'][k]['S_relief']:.3f} (rej {row['arms'][k]['n_reject_relief']})" for k in ("1e4", "1e6")), flush=True)
    pdbs = [r["pdb"] for r in rows]; folds = L.folds_of(pdbs)
    anchor = np.array([r["anchor"] for r in rows])
    out = {"n": len(rows)}
    for T in PR.THRESHOLDS:
        k = PR.tkey(T)
        arms = [r["arms"][k] for r in rows]
        res = {"n_reject_relief": L.mean_se([a["n_reject_relief"] for a in arms]),
               "n_reject_raw": L.mean_se([a["n_reject_raw"] for a in arms]),
               "n_S_empty": int(sum(a["S_empty"] for a in arms))}
        for arm in ("S_relief", "S_raw", "RANDS", "PERMS"):
            v = np.array([a[arm] for a in arms])
            c = ST.compare(v, anchor, folds, names=pdbs, label=f"point_cloud [all, n=126] {arm}@{k} minus anchor")
            print(ST.fmt(c)); res[f"{arm}_vs_anchor"] = c
        for arm, ctrl in (("S_relief", "RANDS"), ("S_relief", "PERMS"), ("S_relief", "S_raw")):
            v = np.array([a[arm] for a in arms]); w = np.array([a[ctrl] for a in arms])
            c = ST.compare(v, w, folds, names=pdbs, label=f"point_cloud [all, n=126] {arm}@{k} minus {ctrl}@{k}")
            print(ST.fmt(c)); res[f"{arm}_vs_{ctrl}"] = c
        out[k] = res
    L.save(OUT, {"what": "rotamer_relief B(iii): arm S on the relieved single point, point cloud", "rows": rows,
                 "report": out}, rows=rows, complete_keys=("pdb", "anchor", "arms"), n_expected=126, module_file=__file__)
    return out


if __name__ == "__main__":
    main()
