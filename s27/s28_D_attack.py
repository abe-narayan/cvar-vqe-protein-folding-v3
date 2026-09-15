#!/usr/bin/env python
"""s27/s28_D_attack.py -- lane D: the standard attack kit on any per-target arm table.

Given a JSONL (or JSON with a `rows` list) of per-target rows and the names of the arm column,
the pdb column and the value column, this prints for every arm against a comparator (default:
production, the DIS built chain from `s27/results/chain_rows.jsonl`, or the point-cloud
`rmsd_dis75` from `pool_rows.jsonl`):

  - `ST.compare` verbatim (paired, SE, MDE, iid and fold CI, W/L, concentration null);
  - the FAIL18 / other-108 split, each stratum its own `ST.compare` (a gain carried by 18
    targets is a different claim);
  - for a FAMILY of arms (a regex over arm names), `ST.best_of_k_within` on the (target x arm)
    matrix: split-half transfer, k_eff, share accounted; the per-target minimum is an order
    statistic, never a result.

Usage:
  python s27/s28_D_attack.py --rows s27/results/s28_C_readout_chain_rows.jsonl \
      --arm-col arm --value-col rmsd_chain --basis chain [--family "TRIM\\[CONS,"] [--arms A B]
  python s27/s28_D_attack.py --rows s27/results/s28_B_rows.jsonl --arm-expr "source|seed|graph|J" ...

ORACLE: every value here is an RMSD to a native (the endpoint); nothing is chosen by it.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402

RESULTS = os.path.join(HERE, "results")
FAIL18 = set(I.FAIL18)


def load_rows(path):
    if path.endswith(".jsonl"):
        with open(path, encoding="utf-8") as f:
            return [json.loads(l) for l in f if l.strip()]
    z = json.load(open(path, encoding="utf-8"))
    return z["rows"] if isinstance(z, dict) else z


def production(basis):
    """Per-target production endpoint: built chain (chain_rows DIS) or point cloud (pool_rows)."""
    out = {}
    if basis == "chain":
        for r in load_rows(os.path.join(RESULTS, "chain_rows.jsonl")):
            if r["config"] == "DIS":
                out[r["pdb"]] = float(r["rmsd_chain"])
    else:
        for r in load_rows(os.path.join(RESULTS, "pool_rows.jsonl")):
            if r["config"] == "DIS":
                out[r["pdb"]] = float(r["rmsd"])
    return out


def arm_name(r, arm_col, arm_expr):
    if arm_expr:
        return "|".join(str(r[k]) if not isinstance(r[k], float) else ("%g" % r[k]) for k in arm_expr.split("|"))
    return str(r[arm_col])


def get_value(r, value_col):
    v = r
    for part in value_col.split("."):
        v = v[part]
    return float(v)


def table(rows, arm_col, arm_expr, value_col):
    T = {}
    for r in rows:
        try:
            v = get_value(r, value_col)
        except (KeyError, TypeError):
            continue
        T.setdefault(arm_name(r, arm_col, arm_expr), {})[r["pdb"]] = v
    return T


def attack_one(name, arm, prod, label_prefix=""):
    pdbs = sorted(p for p in arm if p in prod and np.isfinite(arm[p]))
    if len(pdbs) < 20:
        print(f"  {name}: only {len(pdbs)} paired finite targets, skipped")
        return None
    a = np.array([arm[p] for p in pdbs]); b = np.array([prod[p] for p in pdbs])
    folds = ST.pinned_folds(pdbs)
    r = ST.compare(a, b, folds, names=pdbs, label=f"{label_prefix}{name} vs production (n={len(pdbs)})")
    print(ST.fmt(r))
    out = {"all": {k: v for k, v in r.items() if k != "concentration"}, "n": len(pdbs)}
    out["all"]["concentration"] = r["concentration"]
    for tag, sel in (("FAIL18", [p in FAIL18 for p in pdbs]), ("other108", [p not in FAIL18 for p in pdbs])):
        m = np.array(sel)
        if m.sum() < 5:
            continue
        rs = ST.compare(a[m], b[m], folds[m], names=[p for p, s in zip(pdbs, sel) if s],
                        label=f"    {name} vs production, {tag} (n={int(m.sum())})")
        print(ST.fmt(rs))
        out[tag] = {k: v for k, v in rs.items() if k != "concentration"}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", required=True)
    ap.add_argument("--arm-col", default="arm")
    ap.add_argument("--arm-expr", default=None, help="pipe-joined columns forming the arm name")
    ap.add_argument("--value-col", default="rmsd_chain")
    ap.add_argument("--basis", choices=["chain", "cloud"], default="chain")
    ap.add_argument("--prod-arm", default=None, help="use this arm from the rows as the comparator instead of S27 production")
    ap.add_argument("--arms", nargs="*", default=None)
    ap.add_argument("--family", default=None, help="regex over arm names: best_of_k_within on the matrix")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rows = load_rows(a.rows)
    T = table(rows, a.arm_col, a.arm_expr, a.value_col)
    prod = T[a.prod_arm] if a.prod_arm else production(a.basis)
    names = a.arms or sorted(T)
    res = {"rows": a.rows, "value_col": a.value_col, "basis": a.basis, "prod": a.prod_arm or "S27 production", "arms": {}}
    for nm in names:
        if nm not in T or nm == a.prod_arm:
            continue
        r = attack_one(nm, T[nm], prod)
        if r:
            res["arms"][nm] = r
    if a.family:
        fam = [nm for nm in sorted(T) if re.search(a.family, nm)]
        pdbs = sorted(set.intersection(*(set(T[nm]) for nm in fam)) & set(prod))
        M = np.array([[T[nm][p] for nm in fam] for p in pdbs])
        ok = np.isfinite(M).all(1)
        M = M[ok]
        bk = ST.best_of_k_within(M)
        print(f"\n  FAMILY {a.family!r}: {len(fam)} arms x {int(ok.sum())} targets")
        print(f"    per-target min - mean {bk['observed_gain']:+.4f}   across-target null {bk['null_across_targets']:+.4f}"
              f"   share accounted {bk['share_accounted']:.2f}   k_eff {bk['k_eff']:.2f}")
        print(f"    split-half transfer {bk['split_half']:+.4f} ({100*bk['split_half_frac']:.0f}% of oracle)   {bk['verdict']}")
        print("    arms:", ", ".join(f"{nm} ({c:.0f})" for nm, c in zip(fam, bk["argmin_counts"])))
        res["family"] = {"arms": fam, **{k: v for k, v in bk.items()}}
    if a.out:
        ST.save_atomic(a.out, res, module_file=__file__)
        print("wrote", a.out)


if __name__ == "__main__":
    main()
