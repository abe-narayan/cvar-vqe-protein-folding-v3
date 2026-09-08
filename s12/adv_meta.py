"""ADVERSARIAL AUDIT 3 -- the fibril/lasso enrichment (10/18 vs 6/108, Fisher p=1.2e-06).

Three attacks, all run:
  (a) MULTIPLICITY.  The reported statistic is a UNION of two flags chosen out of a family
      of 11 header flags (forensics_part1.py:151) plus 7 diagnostic flags (:154) plus the
      8+7 flags of fail_headers.py:98.  A union of two chosen post hoc is drawn from
      11 singles + 55 pairwise unions = 66 hypotheses.  Corrected by a MAX-STATISTIC
      PERMUTATION over that whole family (labels shuffled, family minimum p recorded).
  (b) LEAVE-ONE-OUT.  Drop each of the 126 targets in turn; does the enrichment survive?
      Also drop each of the 10 fibril/lasso FAIL18 members one at a time.
  (c) THE LABELS THEMSELVES.  I re-derive fibril/lasso from the raw PDB headers with my
      own reading and diff against forensics'.
"""
from __future__ import annotations
import os, re, sys, json, itertools
import numpy as np
from scipy.stats import fisher_exact

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I

PDBDIRS = [os.path.join(ROOT, "pdbs_ext"), os.path.join(ROOT, "pdbs")]

# my own, independent readings
RE = {
    "fibril": re.compile(r"AMYLOID|FIBRIL|STERIC ZIPPER|PRION|PROTOFIB", re.I),
    "lasso": re.compile(r"LASSO", re.I),
    "cyclic": re.compile(r"CYCLIC|CYCLOTIDE|MACROCYCL|HEAD-TO-TAIL", re.I),
    "membrane": re.compile(r"MICELLE|BICELLE|\bSDS\b|\bDPC\b|MEMBRANE|LIPID|VESICL|TRANSMEMBRANE", re.I),
    "cosolvent": re.compile(r"\bTFE\b|TRIFLUOROETHANOL|DMSO|METHANOL|ACETONITRILE", re.I),
    "bound": re.compile(r"\bBOUND (TO|WITH)\b|\bIN COMPLEX WITH\b|COMPLEXED|ANTAGONIST|INHIBITOR", re.I),
    "designed": re.compile(r"DESIGN|DE NOVO|SYNTHETIC CONSTRUCT", re.I),
    "xray": re.compile(r"X-RAY|ELECTRON|NEUTRON", re.I),
}


def read(pdb):
    p = None
    for d in PDBDIRS:
        q = os.path.join(d, f"{pdb}.pdb")
        if os.path.exists(q):
            p = q
            break
    if p is None:
        return None
    with open(p, "r", errors="ignore") as fh:
        lines = fh.readlines()
    tag = lambda t: "\n".join(l.rstrip() for l in lines if l.startswith(t))
    text = "\n".join([tag("HEADER"), tag("TITLE"), tag("COMPND"), tag("KEYWDS"), tag("EXPDTA")])
    rem = "\n".join(l.rstrip() for l in lines if l.startswith("REMARK 210")
                    or l.startswith("REMARK 245") or l.startswith("REMARK 280"))
    d = {"pdb": pdb, "text": text, "rem": rem,
         "expdta": " ".join(l[10:].strip() for l in lines if l.startswith("EXPDTA")),
         "title": " ".join(l[10:].strip() for l in lines if l.startswith("TITLE")),
         "n_ssbond": sum(1 for l in lines if l.startswith("SSBOND")),
         "n_link": sum(1 for l in lines if l.startswith("LINK"))}
    ctx = text + "\n" + rem
    for k, r in RE.items():
        src = d["expdta"] if k == "xray" else (ctx if k in ("membrane", "cosolvent") else text)
        d[k] = bool(r.search(src))
    d["covalent"] = bool(d["n_ssbond"] or d["n_link"])
    d["nonaqueous"] = d["cosolvent"]
    d["any_flag"] = any(d[k] for k in ("fibril", "lasso", "cyclic", "membrane",
                                       "cosolvent", "bound", "designed", "covalent"))
    return d


def fisher_p(flag, fail):
    a = int((flag & fail).sum()); b = int((flag & ~fail).sum())
    return float(fisher_exact([[a, int(fail.sum()) - a],
                               [b, int((~fail).sum()) - b]])[1]), a, b


def main(n_perm=20000, seed=0):
    tg = I.targets()
    recs = [read(t["pdb"]) for t in tg]
    miss = [t["pdb"] for t, r in zip(tg, recs) if r is None]
    fail = np.array([t["pdb"] in I.FAIL18 for t in tg], bool)
    names = [t["pdb"] for t in tg]

    SINGLES = ["fibril", "lasso", "cyclic", "membrane", "cosolvent", "bound", "designed",
               "xray", "covalent", "nonaqueous", "any_flag"]
    F = {k: np.array([bool(r[k]) if r else False for r in recs]) for k in SINGLES}

    # ---- (c) cross-check against forensics' own labels
    with open(os.path.join(I.RESULTS, "fail_headers.json")) as fh:
        fh_json = json.load(fh)["per_target"]
    diffs = []
    for t in tg:
        o = fh_json.get(t["pdb"], {})
        for k in ("fibril", "lasso", "cyclic", "membrane", "bound"):
            if k in o and bool(o[k]) != bool(F[k][names.index(t["pdb"])]):
                diffs.append((t["pdb"], k, bool(o[k]), bool(F[k][names.index(t["pdb"])])))

    # ---- the family the composite was drawn from: 11 singles + all 55 pairwise unions
    family = {k: F[k] for k in SINGLES}
    for a, b in itertools.combinations(SINGLES, 2):
        family[f"{a}|{b}"] = F[a] | F[b]
    ps = {k: fisher_p(v, fail) for k, v in family.items()}
    reported = ps["fibril|lasso"]

    # ---- (a) max-statistic permutation over the family
    rng = np.random.default_rng(seed)
    keys = list(family)
    M = np.array([family[k] for k in keys])            # (nfam, 126)
    nf = int(fail.sum()); N = len(fail)
    minp = np.empty(n_perm)
    for b in range(n_perm):
        perm = np.zeros(N, bool); perm[rng.choice(N, nf, replace=False)] = True
        best = 1.0
        for row in M:
            a11 = int((row & perm).sum()); a10 = int((row & ~perm).sum())
            p = fisher_exact([[a11, nf - a11], [a10, (N - nf) - a10]])[1]
            if p < best:
                best = p
        minp[b] = best
    p_fwer = float((minp <= reported[0]).mean())
    # also the honest correction if only the 11 singles had been examined
    minp11 = None

    # ---- (b) leave-one-out
    lf = F["fibril"] | F["lasso"]
    loo = []
    for k in range(N):
        m = np.ones(N, bool); m[k] = False
        p, a, bb = fisher_p(lf[m], fail[m])
        loo.append((names[k], p, a, bb))
    loo_max = max(loo, key=lambda x: x[1])
    # drop each fibril/lasso FAIL18 member
    drop_fl = [(names[k], fisher_p(lf[np.arange(N) != k], fail[np.arange(N) != k])[0])
               for k in range(N) if lf[k] and fail[k]]

    out = dict(
        missing_pdbs=miss,
        label_diffs_vs_forensics=diffs,
        reported=dict(p=reported[0], n_fail_flagged=reported[1], n_other_flagged=reported[2]),
        family_size=len(family),
        family_top=sorted(((k, v[0], v[1], v[2]) for k, v in ps.items()), key=lambda x: x[1])[:12],
        bonferroni_p_over_family=float(min(1.0, reported[0] * len(family))),
        bonferroni_p_over_11_singles=float(min(1.0, reported[0] * len(SINGLES))),
        maxstat_permutation_p=p_fwer, n_perm=n_perm,
        loo_worst=dict(dropped=loo_max[0], p=loo_max[1]),
        loo_all_p_max=float(max(x[1] for x in loo)),
        drop_each_flagged_fail=drop_fl,
        singles_p={k: ps[k][0] for k in SINGLES},
        singles_counts={k: [ps[k][1], ps[k][2]] for k in SINGLES},
    )
    print(json.dumps({k: v for k, v in out.items() if k != "family_top"}, indent=1))
    print("family top 12:", json.dumps(out["family_top"], indent=1))
    I.write("adv_meta", out)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 20000)
