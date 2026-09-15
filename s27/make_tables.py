#!/usr/bin/env python
"""s27/make_tables.py -- render every S27 result table as markdown from the summary JSONs
(`pool_summary.json`, `wave2_summary.json`, `wave3_summary.json`, `arms_summary.json`,
`trainability.json`, `redundancy.json`) into `s27/results/tables.md`.  No number is typed by
hand: the report cites these tables.
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results")


def L(name):
    with open(os.path.join(R, name), encoding="utf-8") as fh:
        return json.load(fh)


def ci(c):
    return f"[{c[0]:+.3f}, {c[1]:+.3f}]" if c else "n/a"


def verdict_word(v, x):
    """Map the library verdict + effect to the S27 vocabulary."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "reference"
    if v.startswith("BETTER"):
        return "PROMISING"
    if v.startswith("WORSE"):
        return "WORSE"
    if abs(x) < 0.7:
        return "NULL"
    return "TYPE-M ZONE"


def main():
    pool = L("pool_summary.json"); w2 = L("wave2_summary.json"); w3 = L("wave3_summary.json")
    red = L("redundancy.json"); tr = L("trainability.json")
    arms = L("arms_summary.json") if os.path.exists(os.path.join(R, "arms_summary.json")) else None
    out = []
    C = pool["configs"]
    fam = json.load(open(os.path.join(HERE, "family.json"))) if os.path.exists(os.path.join(HERE, "family.json")) else {}

    out.append(f"## T1. Every configuration, ranked by mean point-cloud RMSD of the top-75 average (n = {pool['n']}; DIS = {pool['dis_mean']:.4f}, random-75 = {pool['rand_mean']:.4f})\n")
    out.append("| rank | configuration | kind | mean RMSD (A) | vs DIS (A) | x MDE | fold-clustered 95% CI | W/L | vs random-75 (A) | verdict |")
    out.append("|---:|---|---|---:|---:|---:|---|---|---:|---|")
    for i, r in enumerate(pool["ranked"]):
        v = C[r["config"]]; a = v["vs_dis"]; b = v["vs_random"]
        out.append(f"| {i+1} | `{r['config']}` | {v['kind']} | {v['mean']:.4f} | {a['effect']:+.4f} | {a['x_mde'] if a['x_mde'] is not None else float('nan'):+.2f} | {ci(a['ci_fold'])} | {a['W']}/{a['L']} | {b['effect']:+.4f} | {verdict_word(a['verdict'], a['x_mde'])} |")
    out.append("")
    out.append("## T2. The singles and composites with their ORACLE diagnostics (rho = Spearman with the candidate RMSD, whole pool and in-band < 3 A; rho_DIS = rank correlation with the shipped score; overlap = share of DIS's top-75 in the channel's top-75; cost per evaluation)\n")
    out.append("| channel | family | mean RMSD | vs DIS | x MDE | vs random | x MDE | rho pool | rho in-band | rho_DIS | overlap | tie frac | ms/eval |")
    out.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    from s27 import ham_lib as HL
    for k, v in sorted([(k, v) for k, v in C.items() if v["kind"] in ("single", "composite")], key=lambda kv: kv[1]["mean"]):
        a = v["vs_dis"]; b = v["vs_random"]
        f = HL.FAMILY.get(k, "legacy term" if k.startswith("LEG_") else "composite")
        out.append(f"| `{k}` | {f} | {v['mean']:.4f} | {a['effect']:+.4f} | {a['x_mde'] if a['x_mde'] is not None else 0:+.2f} | {b['effect']:+.4f} | {b['x_mde'] if b['x_mde'] is not None else 0:+.2f} | {v['rho_pool']:+.3f} | {v['rho_inband']:+.3f} | {v.get('rho_dis', float('nan')):+.3f} | {v['overlap_dis']:.2f} | {v['tie_frac_top']:.2f} | {v.get('cost_ms', float('nan')):.3f} |")
    out.append("")
    out.append("## T3. Equal-weight pairs DIS + X: against DIS and against the rank-permuted control (X's marginal kept, its correspondence destroyed)\n")
    out.append("| pair | mean | vs DIS | x MDE | fold CI | vs permuted | x MDE | fold CI (perm) | reading |")
    out.append("|---|---:|---:|---:|---|---:|---:|---|---|")
    for k, v in sorted([(k, v) for k, v in C.items() if v["kind"] == "pair" and "vs_perm" in v], key=lambda kv: kv[1]["vs_perm"]["effect"]):
        a = v["vs_dis"]; p = v["vs_perm"]
        rd = ("carries information (below its noise)" if p["ci_fold"][1] < 0 else "worse than its own noise" if p["ci_fold"][0] > 0 else "indistinguishable from noise")
        out.append(f"| `{k}` | {v['mean']:.4f} | {a['effect']:+.4f} | {a['x_mde'] or 0:+.2f} | {ci(a['ci_fold'])} | {p['effect']:+.4f} | {p['x_mde'] or 0:+.2f} | {ci(p['ci_fold'])} | {rd} |")
    out.append("")
    out.append("## T4. Regulariser weights (H3): the weight chosen leave-fold-out, the held-out effect, and the grid priced as an order statistic\n")
    out.append("| complement | held-out effect | x MDE | fold CI | weights chosen per fold | full-leakage best | k_eff |")
    out.append("|---|---:|---:|---|---|---:|---:|")
    for k, v in sorted(pool["nested_weight"].items(), key=lambda kv: kv[1]["held_out"]["effect"]):
        h = v["held_out"]
        out.append(f"| `{k}` | {h['effect']:+.4f} | {h['x_mde'] or 0:+.2f} | {ci(h['ci_fold'])} | {list(v['chosen_w_per_fold'].values())} | {v['full_leak_best']:+.4f} | {v['best_of_k_within'].get('k_eff', float('nan')):.2f} |")
    n4 = w3["nested_SS_MATCH_4weights"]
    out.append(f"\nSS_MATCH with the four-weight grid (0.25, 0.5, 0.75, 1.0): held-out {n4['held_out_effect']:+.4f} A ({n4['x_mde']:+.2f}x MDE), fold CI {ci(n4['ci_fold'])}, chosen {n4['chosen']}, `best_of_k_within`: {n4['best_of_k'].get('verdict')} (k_eff {n4['best_of_k'].get('k_eff'):.2f}).\n")
    out.append("## T5. The m-ladder (H9): mean point-cloud RMSD of the top-m average; in parentheses the paired effect against DIS at the same m (a value beyond its MDE with the fold CI excluding zero is marked *)\n")
    MS = [3, 5, 10, 25, 50, 75, 100, 150]
    arms_l = ["DIS", "DISTPOT", "CONS", "DIS+DISTPOT", "DIS+CONS", "DIS+ENV", "DIS+SS_MATCH", "DIS+0.5*SS_MATCH"]
    out.append("| m | random-m | " + " | ".join(f"`{a}`" for a in arms_l) + " |")
    out.append("|---:|---:|" + "---:|" * len(arms_l))
    for m in MS:
        cells = []
        for a in arms_l:
            r = w2[f"H9|{a}|{m}"]
            star = "*" if (r["ci_fold"][1] < 0 and abs(r["x_mde"]) >= 1) else ""
            cells.append(f"{r['mean']:.3f} ({r['effect']:+.3f}){star}")
        out.append(f"| {m} | {w2[f'H9|DIS|{m}']['rand_m']:.3f} | " + " | ".join(cells) + " |")
    out.append("")
    out.append("## T6. In-band re-ranking (H10): DIS's top-150, then 75 kept by X\n")
    out.append("| rule | mean | vs DIS top-75 | x MDE | fold CI | W/L |")
    out.append("|---|---:|---:|---:|---|---|")
    for k, v in sorted([(k, v) for k, v in w2.items() if k.startswith("H10|")], key=lambda kv: kv[1]["mean"]):
        out.append(f"| `{k[4:]}` | {v['mean']:.4f} | {v['effect']:+.4f} | {v['x_mde']:+.2f} | {ci(v['ci_fold'])} | {v['W']}/{v['L']} |")
    out.append("")
    out.append("## T7. Non-additive forms (H11), near-corpus potentials (H12) and the SS_MATCH variants (H13)\n")
    out.append("| configuration | hypothesis | mean | vs DIS | x MDE | fold CI | folds same sign | W/L |")
    out.append("|---|---|---:|---:|---:|---|---|---|")
    for k, v in sorted([(k, v) for k, v in w2.items() if k.startswith("H11|") or k.startswith("H12|")], key=lambda kv: kv[1]["mean"]):
        out.append(f"| `{k.split('|',1)[1]}` | {k.split('|')[0]} | {v['mean']:.4f} | {v['effect']:+.4f} | {v['x_mde']:+.2f} | {ci(v['ci_fold'])} | {v['folds_same_sign']}/5 | {v['W']}/{v['L']} |")
    for k, v in sorted([(k, v) for k, v in w3.items() if k != "nested_SS_MATCH_4weights" and k != "DIS"], key=lambda kv: kv[1]["mean"]):
        out.append(f"| `{k}` | H13 | {v['mean']:.4f} | {v['effect']:+.4f} | {v['x_mde']:+.2f} | {ci(v['ci_fold'])} | {v['folds_same_sign']}/5 | {v['W']}/{v['L']} |")
    out.append("")
    out.append("## T8. Complementarity and redundancy: partial Spearman of each channel with the ORACLE candidate RMSD given DIS (positive = adds correct ranking information), and its rank correlation with DIS\n")
    out.append("| channel | partial rho given DIS | rho with DIS | DIS+X effect on the top-75 average |")
    out.append("|---|---:|---:|---:|")
    for k, p in sorted(red["partial_rho_given_dis"].items(), key=lambda kv: -kv[1]):
        e = C.get(f"DIS+{k}", {}).get("vs_dis", {}).get("effect")
        out.append(f"| `{k}` | {p:+.3f} | {red['rho_with_dis'][k]:+.3f} | {'' if e is None else f'{e:+.4f}'} |")
    out.append("")
    if arms:
        out.append("## T9. The genuine CVaR-VQE arm (9 qubits, alpha 0.18, T 0.5, 80 Adam iterations, exact parameter-shift gradient), point cloud: each configuration's VQE selection against DIS's VQE selection, seed 0, with the seed-1 replication, the realised tail m, the set-equality gate, and eps = RMSD_VQE - RMSD_top-m\n")
        out.append("| configuration | VQE s0 | vs DIS-VQE | x MDE | fold CI | VQE s1 | effect s1 | s1 inside s0's CI | m | max eps | entropy bits | ESS | gate |")
        out.append("|---|---:|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---|")
        for cfg, ds in sorted(arms["vqe"].items(), key=lambda kv: kv[1]["0"]["mean_vqe"]):
            d0 = ds["0"]; a = d0["vs_dis_vqe"]; d1 = ds.get("1")
            out.append(f"| `{cfg}` | {d0['mean_vqe']:.4f} | {a['effect']:+.4f} | {a['x_mde'] or 0:+.2f} | {ci(a['ci_fold'])} | {d1['mean_vqe'] if d1 else float('nan'):.4f} | {d1['vs_dis_vqe']['effect'] if d1 else float('nan'):+.4f} | {'yes' if ds.get('replicates') else 'no'} | {d0['mean_m']:.1f} | {d0['eps_max']:.1e} | {d0['entropy_bits']:.2f} | {d0['ess']:.1f} | {'pass' if d0['gate_pass_all'] else 'FAIL'} |")
        out.append("")
        if "chain" in arms:
            out.append("## T10. The built-chain endpoint (the production projection of the same top-75 averages) against DIS, with the same sets on the point cloud beside it\n")
            out.append("| configuration | built chain | vs DIS | x MDE | fold CI | W/L | point cloud | vs DIS | x MDE |")
            out.append("|---|---:|---:|---:|---|---|---:|---:|---:|")
            for cfg, d in sorted(arms["chain"].items(), key=lambda kv: kv[1]["mean_chain"]):
                a = d["vs_dis_chain"]; b = d["vs_dis_cloud"]
                out.append(f"| `{cfg}` | {d['mean_chain']:.4f} | {a['effect']:+.4f} | {a['x_mde'] or 0:+.2f} | {ci(a['ci_fold'])} | {a['W']}/{a['L']} | {d['mean_cloud']:.4f} | {b['effect']:+.4f} | {b['x_mde'] or 0:+.2f} |")
            out.append("")
    out.append("## T11. Trainability (H8): gradient variance of the deployed objective at n = 9, depth 3, alpha 0.18, T 0.5, 120 theta draws, median over 12 targets, under the deployed rank currency and under a bounded standardisation that keeps the channel's own gaps\n")
    out.append("| channel | Var[dF/dtheta_0], zrank | Var, asinh-MAD | ratio to DIS (asinh) |")
    out.append("|---|---:|---:|---:|")
    for c in tr["summary"]:
        if c.endswith("|zrank"):
            n = c.split("|")[0]
            out.append(f"| `{n}` | {tr['summary'][c]['var_g0_median']:.4e} | {tr['summary'][n+'|asinh']['var_g0_median']:.4e} | {tr['summary'][n+'|asinh']['ratio_to_DIS_same_std']:.2f} |")
    out.append("")
    with open(os.path.join(R, "tables.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))
    print("wrote", os.path.join(R, "tables.md"), len(out), "lines")


if __name__ == "__main__":
    main()
