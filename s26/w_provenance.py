#!/usr/bin/env python
"""s26/w_provenance.py -- window provenance: whole-peptide / terminal / interior / fragment
windows in the K = 500 pool and the shipped top-75.  Lane W, Sprint 26.  Pre-registration:
s26/PREREG_window_provenance.md.

`census` is native-free (codes and sequences only; the universe is loaded blind).  `oracle`
(GATED) reads the universe's `rr` per member (ORACLE RMSD-to-native of single windows) and
contrasts the classes inside the top-75 and the pool.  `readout` (GATED, only if the ORACLE
class contrast inside the top-75 clears its MDE) re-weights the shipped top-75 by class.

    python s26/jobrun.py --agent W --tag CPU --name w_provenance_census --est-ram 0.3 -- python s26/w_provenance.py census
    python s26/jobrun.py --agent W --tag CPU --name w_provenance_oracle --est-ram 0.3 -- python s26/w_provenance.py oracle
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (ROOT, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s12 import instrument as I                      # noqa: E402
from s24 import stats_lib as ST                      # noqa: E402
from core import data as D                           # noqa: E402
import w_selfcopy as W                               # noqa: E402

K, TOPM = W.K, W.TOPM
CLASSES = ("whole", "terminal", "interior", "fragment", "unresolved")
_IDX = {}


def substring_index(n):
    """length-n window string -> (class, parent kind) over the 787 peptides then the 6,003
    fragments; the FIRST parent in database order wins, and multi-parent windows are counted."""
    if n in _IDX:
        return _IDX[n]
    idx, multi = {}, 0
    for kind, chains in (("peptide", D.load()), ("fragment", D.load_fragments())):
        for p in chains:
            s = p.seq; m = len(s)
            for off in range(0, m - n + 1):
                w = s[off:off + n]
                if kind == "peptide":
                    cls = "whole" if m == n else ("terminal" if off in (0, m - n) else "interior")
                else:
                    cls = "fragment"
                if w in idx:
                    multi += 1
                else:
                    idx[w] = (cls, kind)
    _IDX[n] = (idx, multi)
    return _IDX[n]


def classify(u, seq, members):
    idx, _ = substring_index(len(seq))
    S = np.asarray(u["S"]); org = np.asarray(u["org"], bool)
    out = []
    for k in members:
        w = "".join(D.ALPHABET[int(c)] for c in S[int(k)])
        hit = idx.get(w)
        if hit is None:
            out.append("unresolved")
        elif org[int(k)] and hit[1] == "peptide":
            out.append(hit[0])
        elif (not org[int(k)]) and hit[1] == "fragment":
            out.append("fragment")
        else:
            # the string exists in the other bank as well; trust the universe's org flag
            out.append("fragment" if not org[int(k)] else hit[0] if hit[1] == "peptide" else "unresolved")
    return out


def census(verbose=True):
    tg = I.targets()
    rows = []
    for t in tg:
        pdb, seq, n = t["pdb"], t["seq"], t["n"]
        u = W.load_blind(pdb)
        order = np.asarray(u["order"], int); pool = order[:K]
        rec = I.shipped_record(pdb); sub = np.asarray(rec["sub"], int)
        cls_pool = classify(u, seq, pool)
        cls_top = [cls_pool[int(p)] for p in sub]
        row = {"pdb": pdb, "n": n, "fold": t["fold"],
               "pool": {c: int(sum(1 for x in cls_pool if x == c)) for c in CLASSES},
               "top75": {c: int(sum(1 for x in cls_top if x == c)) for c in CLASSES},
               "class_of_pool_member": cls_pool}
        rows.append(row)
        if verbose:
            print("  %s n=%2d pool %s | top75 %s" % (pdb, n, {c: row["pool"][c] for c in CLASSES if row["pool"][c]},
                                                      {c: row["top75"][c] for c in CLASSES if row["top75"][c]}), flush=True)
        del u
    summ = {"n": len(rows), "multi_parent_windows_by_length": {n: substring_index(n)[1] for n in sorted({r["n"] for r in rows})},
            "pool_frac_mean": {c: float(np.mean([r["pool"][c] / K for r in rows])) for c in CLASSES},
            "top75_frac_mean": {c: float(np.mean([r["top75"][c] / TOPM for r in rows])) for c in CLASSES},
            "targets_with_whole_in_top75": int(sum(1 for r in rows if r["top75"]["whole"])),
            "targets_with_terminal_in_top75": int(sum(1 for r in rows if r["top75"]["terminal"])),
            "unresolved_total": int(sum(r["pool"]["unresolved"] for r in rows))}
    p = W.save("provenance_census", {"label": "window provenance census, native-free", "summary": summ, "rows": rows},
               rows=rows, n_expected=len(tg), complete_keys=("pdb", "pool", "top75", "class_of_pool_member"))
    print(json.dumps(summ, indent=1)); print("  ->", p)
    return rows


def oracle(verbose=True):
    W.require_signoff("provenance oracle")
    Z = W.load_result("provenance_census")
    if not (Z and Z.get("complete")):
        raise SystemExit("run `census` first")
    rows = []
    for r in Z["rows"]:
        u = I.load_univ(r["pdb"]); order = np.asarray(u["order"], int); pool = order[:K]
        rr = np.asarray(u["rr"], float)[pool]
        sub = np.asarray(I.shipped_record(r["pdb"])["sub"], int)
        cls = np.array(r["class_of_pool_member"])
        pep = np.isin(cls, ["whole", "terminal", "interior"])
        row = {"pdb": r["pdb"], "fold": r["fold"], "pool": {}, "top75": {}}
        for c in CLASSES:
            m = cls == c
            row["pool"][c] = {"n": int(m.sum()), "mean_rr": float(rr[m].mean()) if m.any() else None,
                              "best_rr": float(rr[m].min()) if m.any() else None}
            mt = m[sub]
            row["top75"][c] = {"n": int(mt.sum()), "mean_rr": float(rr[sub][mt].mean()) if mt.any() else None}
        for name, m in (("peptide_any", pep), ("whole_or_terminal", np.isin(cls, ["whole", "terminal"]))):
            row["pool"][name] = {"n": int(m.sum()), "mean_rr": float(rr[m].mean()) if m.any() else None}
            mt = m[sub]
            row["top75"][name] = {"n": int(mt.sum()), "mean_rr": float(rr[sub][mt].mean()) if mt.any() else None}
        rows.append(row)
        del u
    out = {"label": "ORACLE class contrast (single-window basis, rr)", "rows": rows, "stats": {}}
    for where in ("top75", "pool"):
        for a, b in (("whole_or_terminal", "fragment"), ("peptide_any", "fragment"), ("interior", "fragment"), ("whole_or_terminal", "interior")):
            sel = [r for r in rows if r[where][a]["mean_rr"] is not None and r[where][b]["mean_rr"] is not None]
            if len(sel) < 10:
                continue
            pd = [r["pdb"] for r in sel]; folds = ST.pinned_folds(pd)
            res = ST.compare(np.array([r[where][a]["mean_rr"] for r in sel]), np.array([r[where][b]["mean_rr"] for r in sel]), folds, names=pd,
                             label="%s: class %s minus class %s, per-target mean ORACLE rr (single-window basis), n=%d targets with both" % (where, a, b, len(sel)))
            out["stats"]["%s_%s_minus_%s" % (where, a, b)] = res
            if verbose:
                print(ST.fmt(res))
    p = W.save("provenance_oracle", out)
    print("  ->", p)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("census", "oracle"))
    a = ap.parse_args(argv)
    if a.cmd == "census":
        census()
    else:
        oracle()
    return 0


if __name__ == "__main__":
    sys.exit(main())
