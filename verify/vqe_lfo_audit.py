"""Is `VQE_LFO` genuinely leave-fold-out, or was it read off the instrument?

`core.pipeline.VQE_LFO` maps fold -> (alpha, T) and the docstring claims fold f's cell was
selected using ONLY the other four folds -- which is what makes consulting it not tuning
on the 126-target instrument. That claim is checkable: `s8/integrate_vqe.json` holds the
per-target selected RMSD of every (alpha, T) arm together with each target's fold, so the
leave-fold-out choice can be RECOMPUTED here and compared with the shipped table.

Also reported, because it is the number that decides whether the component earns its
place: the LFO table's own mean, the single best global cell (which is what tuning on the
instrument WOULD have picked), and the shipped argmin baseline.
"""
from __future__ import annotations

import json
import os
import re
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

SRC = os.path.join(_ROOT, "s8", "integrate_vqe.json")


def main():
    d = json.load(open(SRC))
    folds = np.asarray(d["folds"], int)
    per = d["per_target"]
    base = float(d["baseline_sel"])
    arms = {k: np.asarray(v, float) for k, v in per.items()
            if k.startswith("vqe_a")}
    out = {"source": os.path.relpath(SRC, _ROOT), "n_targets": d["n_targets"],
           "iters": d["iters"], "n_qubits": d["n_qubits"], "layers": d["layers"],
           "dim": d["dim"], "baseline_sel_shipped_argmin": base,
           "n_arms": len(arms)}

    def parse(name):
        m = re.match(r"vqe_a([0-9.]+)_T([0-9.]+)$", name)
        return (float(m.group(1)), float(m.group(2))) if m else None

    # the single best GLOBAL cell -- what tuning on the instrument would have chosen
    gbest = min(arms, key=lambda k: arms[k].mean())
    out["best_global_cell"] = {"cell": gbest, "params": parse(gbest),
                               "mean": float(arms[gbest].mean())}

    # recompute the leave-fold-out choice: for each fold f, pick the cell with the best
    # mean over targets NOT in f, then score it ON f.
    ufolds = sorted(set(folds.tolist()))
    recomputed, lfo_vals = {}, np.zeros(len(folds))
    for f in ufolds:
        out_mask = folds != f
        in_mask = folds == f
        pick = min(arms, key=lambda k: arms[k][out_mask].mean())
        recomputed[f] = parse(pick)
        lfo_vals[in_mask] = arms[pick][in_mask]
    out["recomputed_LFO"] = {int(k): list(v) for k, v in recomputed.items()}
    out["recomputed_LFO_mean"] = float(lfo_vals.mean())

    from core.pipeline import VQE_LFO
    shipped = {int(k): tuple(float(x) for x in v) for k, v in VQE_LFO.items()}
    out["shipped_VQE_LFO"] = {k: list(v) for k, v in shipped.items()}
    match = {f: (tuple(recomputed[f]) == shipped.get(f)) for f in ufolds}
    out["per_fold_match"] = {int(k): bool(v) for k, v in match.items()}
    out["LFO_TABLE_REPRODUCES"] = all(match.values())

    # score the SHIPPED table the same way, so its mean is comparable
    sh_vals = np.zeros(len(folds))
    ok = True
    for f in ufolds:
        a, T = shipped.get(f, (None, None))
        # arm names are written with the repr the experiment used, so try both forms
        cands = [f"vqe_a{a}_T{T}", f"vqe_a{a:g}_T{T:g}"]
        key = next((c for c in cands if c in arms), None)
        if key is None:
            ok = False
            out["shipped_LFO_key_miss"] = {"fold": int(f), "tried": cands,
                                           "available": sorted(arms)[:8]}
            break
        sh_vals[folds == f] = arms[key][folds == f]
    out["shipped_LFO_mean"] = float(sh_vals.mean()) if ok else None
    out["shipped_LFO_vs_argmin"] = (float(sh_vals.mean() - base) if ok else None)
    out["best_global_vs_argmin"] = float(arms[gbest].mean() - base)
    # the honest gap: LFO must be WORSE than the global best, or it was not held out
    out["LFO_is_not_better_than_global_best"] = (
        None if not ok else bool(sh_vals.mean() >= arms[gbest].mean() - 1e-12))

    print(json.dumps(out, indent=2, sort_keys=True, default=float))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "vqe_lfo_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=float)
    return out

if __name__ == "__main__":
    main()
