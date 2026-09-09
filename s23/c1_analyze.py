"""s23/c1_analyze.py -- aggregate `c1_probweight_raw.json` into the pre-registered comparisons.

Per target: average the 4 seeds (matched seed-for-seed between trained/untrained, per
`c1_probweight.py`'s construction) into one point per (target, T, arm). Then, across the 126
targets, the four paired comparisons the PREREG calls for:

    PRIMARY   W_trained  - C_trained    (weighted vs size-matched classical bar, trained tail)
    CONTROL   W_untrained - C_untrained (same question, on the untrained circuit's own tail)
    MECH      W_trained  - W_untrained  (does training move the WEIGHTED answer)
    MECH2     C_trained  - C_untrained  (does training move the answer even under UNIFORM
                                          weights, i.e. purely via which m gets realised --
                                          s22 D8's mechanism, replicated as a check here)
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

from s23 import qc_lib as QL   # noqa: E402


def load_raw():
    return json.load(open(os.path.join(QL.RESULTS, "c1_probweight_raw.json")))


def per_target_frame(raw, T: str):
    pdbs, folds = [], []
    w_tr, c_tr, w_un, c_un = [], [], [], []
    m_tr, m_un = [], []
    gate_tr, gate_un = [], []
    ent_tr, ent_un = [], []
    for row in raw["rows"]:
        seedrows = row["by_T"][T]
        pdbs.append(row["pdb"]); folds.append(row["fold"])
        w_tr.append(np.mean([s["rmsd_w_trained"] for s in seedrows]))
        c_tr.append(np.mean([s["rmsd_c_trained"] for s in seedrows]))
        w_un.append(np.mean([s["rmsd_w_untrained"] for s in seedrows]))
        c_un.append(np.mean([s["rmsd_c_untrained"] for s in seedrows]))
        m_tr.append(np.mean([s["m_trained"] for s in seedrows]))
        m_un.append(np.mean([s["m_untrained"] for s in seedrows]))
        gate_tr.append(np.mean([1.0 if s["gate1_trained"] else 0.0 for s in seedrows]))
        gate_un.append(np.mean([1.0 if s["gate1_untrained"] else 0.0 for s in seedrows]))
        ent_tr.append(np.mean([s["entropy_bits_trained"] for s in seedrows]))
        ent_un.append(np.mean([s["entropy_bits_untrained"] for s in seedrows]))
    return dict(pdb=np.array(pdbs), fold=np.array(folds, int),
               w_tr=np.array(w_tr), c_tr=np.array(c_tr),
               w_un=np.array(w_un), c_un=np.array(c_un),
               m_tr=np.array(m_tr), m_un=np.array(m_un),
               gate_tr=np.array(gate_tr), gate_un=np.array(gate_un),
               ent_tr=np.array(ent_tr), ent_un=np.array(ent_un))


def analyze_T(raw, T: str) -> dict:
    fr = per_target_frame(raw, T)
    n = len(fr["pdb"])
    out = {"T": T, "n": n,
          "gate1_pass_rate_trained": float(fr["gate_tr"].mean()),
          "gate1_pass_rate_untrained": float(fr["gate_un"].mean()),
          "mean_m_trained": float(fr["m_tr"].mean()), "mean_m_untrained": float(fr["m_un"].mean()),
          "mean_entropy_bits_trained": float(fr["ent_tr"].mean()),
          "mean_entropy_bits_untrained": float(fr["ent_un"].mean()),
          "mean_w_trained": float(fr["w_tr"].mean()), "mean_c_trained": float(fr["c_tr"].mean()),
          "mean_w_untrained": float(fr["w_un"].mean()), "mean_c_untrained": float(fr["c_un"].mean())}
    primary = QL.paired_stats(fr["w_tr"], fr["c_tr"], fr["fold"],
                              label="W_trained - C_trained (a: THE BAR)")
    control = QL.paired_stats(fr["w_un"], fr["c_un"], fr["fold"],
                              label="W_untrained - C_untrained (c's own bar)")
    mech = QL.paired_stats(fr["w_tr"], fr["w_un"], fr["fold"],
                           label="W_trained - W_untrained (does training move the WEIGHTED avg)")
    mech2 = QL.paired_stats(fr["c_tr"], fr["c_un"], fr["fold"],
                            label="C_trained - C_untrained (m-realisation mechanism, uniform wts)")
    ab_equal = QL.paired_stats(fr["c_tr"], fr["c_tr"], fr["fold"], label="(a) vs (b) sanity")
    out["primary"] = {**primary, "verdict": QL.verdict(primary)}
    out["control_c_own_bar"] = {**control, "verdict": QL.verdict(control)}
    out["mechanism_training_effect_on_weighted"] = {**mech, "verdict": QL.verdict(mech)}
    out["mechanism_training_effect_on_classical"] = {**mech2, "verdict": QL.verdict(mech2)}
    return out


def run():
    raw = load_raw()
    assert raw["complete"] and raw["n_expected"] == 126, \
        f"raw not complete: complete={raw['complete']} n={len(raw['rows'])}"
    Ts = [str(raw["config"]["T_primary"])] + [str(t) for t in raw["config"]["T_secondary"]]
    out = {"config": raw["config"], "by_T": {}, "complete": True}
    for T in Ts:
        out["by_T"][T] = analyze_T(raw, T)
    QL.save_json("c1_probweight_analysis.json", out)
    return out


def report(out=None):
    if out is None:
        out = json.load(open(os.path.join(QL.RESULTS, "c1_probweight_analysis.json")))
    print("=== H_C1 -- probability-weighted CVaR-VQE readout ===")
    print(f"alpha={out['config']['alpha']}  T_primary={out['config']['T_primary']}  "
          f"T_secondary={out['config']['T_secondary']}  seeds={out['config']['seeds']}")
    for T, d in out["by_T"].items():
        tag = "PRIMARY" if T == str(out["config"]["T_primary"]) else "secondary"
        print(f"\n--- T={T} [{tag}] ---")
        print(f"  gate1 pass: trained {d['gate1_pass_rate_trained']*100:.1f}%  "
              f"untrained {d['gate1_pass_rate_untrained']*100:.1f}%")
        print(f"  mean m: trained {d['mean_m_trained']:.1f}  untrained {d['mean_m_untrained']:.1f}")
        print(f"  mean entropy (bits, max 9): trained {d['mean_entropy_bits_trained']:.2f}  "
              f"untrained {d['mean_entropy_bits_untrained']:.2f}")
        p = d["primary"]
        print(f"  PRIMARY  W_trained={d['mean_w_trained']:.4f}  C_trained(bar)={d['mean_c_trained']:.4f}")
        print(f"    diff={p['mean_diff']:+.4f}  SE={p['se']:.4f}  MDE={p['mde']:.4f}  "
              f"|eff|/MDE={p['abs_over_mde']:.2f}  CI_iid={p['ci95_iid']}  CI_fold={p['ci95_fold']}  "
              f"W/L={p['n_better']}/{p['n_worse']}  VERDICT={p['verdict']}")
        c = d["control_c_own_bar"]
        print(f"  CONTROL  W_untrained={d['mean_w_untrained']:.4f}  C_untrained(bar)={d['mean_c_untrained']:.4f}")
        print(f"    diff={c['mean_diff']:+.4f}  |eff|/MDE={c['abs_over_mde']:.2f}  "
              f"CI_fold={c['ci95_fold']}  VERDICT={c['verdict']}")
        m1 = d["mechanism_training_effect_on_weighted"]
        print(f"  MECH(W_tr-W_un)  diff={m1['mean_diff']:+.4f}  |eff|/MDE={m1['abs_over_mde']:.2f}  "
              f"CI_fold={m1['ci95_fold']}  VERDICT={m1['verdict']}")
        m2 = d["mechanism_training_effect_on_classical"]
        print(f"  MECH(C_tr-C_un)  diff={m2['mean_diff']:+.4f}  |eff|/MDE={m2['abs_over_mde']:.2f}  "
              f"CI_fold={m2['ci95_fold']}  VERDICT={m2['verdict']}")


if __name__ == "__main__":
    o = run()
    report(o)
