"""S15 AUDIT -- the leave-fold-out discipline of the instrument, checked from data.

Two questions the pinned constants rest on and that nothing in `s12.instrument` asserts:

  L1  Is every one of the 126 tuning targets scored by a model of a fold it is NOT in?
  L2  Is the window universe for each target actually free of that target's own peptide,
      and of every peptide in its fold?

    python -m s15.audit_leak
"""
from __future__ import annotations
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT = os.path.join(ROOT, "s15", "results")
os.makedirs(OUT, exist_ok=True)


def main():
    import peptide_db as db
    from s12 import instrument as I
    folds = db.folds(5)
    clusters = db.clusters()
    peps = list(db.load())
    byseq = {p.seq: p for p in peps}
    tg = I.targets()
    rows = []
    bad_self, bad_fold, bad_cluster = [], [], []
    for t in tg:
        u = I.load_univ(t["pdb"])
        seq = t["seq"]; fold = int(t["fold"])
        # the fold the target's own sequence belongs to, recomputed from the pinned file
        f_true = int(folds[seq])
        # every ORG window in the universe: decode its sequence and check origin
        S = np.asarray(u["S"], int)
        org = np.asarray(u["org"], bool)
        wseqs = set()
        for row in S[org][:200000]:
            wseqs.add("".join(I.ALPHABET[c] for c in row))
        n_self = sum(1 for w in wseqs if w in seq or seq in w)
        # peptide-database members whose windows could be present
        contributors = {p.seq for p in peps if folds[p.seq] != fold and p.seq != seq}
        same_fold_in_univ = [w for w in wseqs if w in byseq and folds[w] == fold]
        same_cluster = [w for w in wseqs
                        if w in byseq and clusters.get(w) == clusters.get(seq)]
        row = {"pdb": t["pdb"], "n": t["n"], "fold_in_npz": fold, "fold_recomputed": f_true,
               "fold_matches": fold == f_true,
               "n_distinct_org_windows": len(wseqs),
               "n_windows_substring_of_target": n_self,
               "n_full_members_from_own_fold": len(same_fold_in_univ),
               "n_full_members_from_own_cluster": len(same_cluster),
               "n_eligible_contributors": len(contributors)}
        if fold != f_true:
            bad_fold.append(t["pdb"])
        if same_fold_in_univ:
            bad_cluster.append((t["pdb"], same_fold_in_univ[:3]))
        if n_self:
            bad_self.append((t["pdb"], n_self))
        rows.append(row)
        if len(rows) % 25 == 0:
            print(f"[{len(rows)}/126]", flush=True)
            with open(os.path.join(OUT, "audit_leak.json"), "w") as fh:
                json.dump({"rows": rows, "complete": False}, fh, indent=1)

    agg = {"n_targets": len(rows),
           "fold_label_mismatches": bad_fold,
           "targets_with_own_fold_members_in_universe": bad_cluster[:20],
           "n_targets_with_own_fold_members": len(bad_cluster),
           "targets_whose_universe_contains_a_window_that_is_a_substring_of_the_target":
               bad_self[:20],
           "n_such_targets": len(bad_self),
           "max_substring_windows": max((n for _, n in bad_self), default=0),
           "note": "org=True windows only (peptide database); protein fragments are the "
                   "other 80% of every pool and are held out by fold, not by identity."}
    print(json.dumps(agg, indent=1))
    with open(os.path.join(OUT, "audit_leak.json"), "w") as fh:
        json.dump({"aggregate": agg, "rows": rows, "n_expected": 126,
                   "complete": len(rows) == 126}, fh, indent=1)
    print("wrote", os.path.join(OUT, "audit_leak.json"))


if __name__ == "__main__":
    main()
