"""FAIL18 forensics, step 1: PDB provenance for all 126 tuning targets.

Reads only headers/records of pdbs_ext/<PDB>.pdb (deposited metadata -- this is a
DIAGNOSTIC / label-audit arm, not a deployable signal).  Classifies each target by
experimental method, number of deposited models, covalent topology (SSBOND / LINK /
HETATM), and keyword-derived context (amyloid/fibril, lasso, cyclic, micelle/bicelle/TFE,
partner-bound).  Compares FAIL18 against the other 108.
"""
from __future__ import annotations
import os, re, sys, json, collections
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I

PDBDIRS = [os.path.join(ROOT, "pdbs_ext"), os.path.join(ROOT, "pdbs")]

FIBRIL_RE = re.compile(r"AMYLOID|FIBRIL|STERIC ZIPPER|PRION|PROTOFIB", re.I)
LASSO_RE = re.compile(r"LASSO", re.I)
CYCLIC_RE = re.compile(r"CYCLIC|CYCLOTIDE|MACROCYCL|HEAD-TO-TAIL", re.I)
MEMB_RE = re.compile(r"MICELLE|BICELLE|SDS|DPC|MEMBRANE|LIPID|TFE|TRIFLUOROETHANOL|"
                     r"VESICL|DODECYLPHOSPHOCHOLINE|TRANSMEMBRANE", re.I)
BOUND_RE = re.compile(r"\bBOUND (TO|WITH)\b|\bIN COMPLEX WITH\b|COMPLEXED", re.I)
FRAG_RE = re.compile(r"^COMPND.*FRAGMENT:", re.I)


def path_for(pdb):
    for d in PDBDIRS:
        p = os.path.join(d, f"{pdb}.pdb")
        if os.path.exists(p):
            return p
    return None


def annotate(pdb):
    p = path_for(pdb)
    rec = dict(pdb=pdb, path=p)
    if p is None:
        return rec
    with open(p, "r", errors="ignore") as fh:
        lines = fh.readlines()
    head = [l.rstrip("\n") for l in lines
            if l[:6].strip() in ("HEADER", "TITLE", "COMPND", "KEYWDS", "EXPDTA",
                                 "NUMMDL", "REMARK", "SSBOND", "LINK", "SEQRES",
                                 "MODRES", "HET")]
    text = "\n".join(l for l in head if l[:6].strip() in
                     ("HEADER", "TITLE", "COMPND", "KEYWDS", "EXPDTA"))
    # remark 210 / 350 sample conditions
    rem = "\n".join(l for l in head if l.startswith("REMARK"))
    expdta = " ".join(l[10:].strip() for l in head if l.startswith("EXPDTA"))
    nummdl = [l for l in head if l.startswith("NUMMDL")]
    nmod = int(nummdl[0][10:].strip()) if nummdl else 1
    title = " ".join(l[10:].strip() for l in head if l.startswith("TITLE"))
    hdr = ([l[10:50].strip() for l in head if l.startswith("HEADER")] or [""])[0]
    kw = " ".join(l[10:].strip() for l in head if l.startswith("KEYWDS"))
    ssb = [l for l in lines if l.startswith("SSBOND")]
    link = [l for l in lines if l.startswith("LINK")]
    het = sorted({l[17:20].strip() for l in lines if l.startswith("HETATM")}
                 - {"HOH", "DOD"})
    modres = sorted({l[12:15].strip() for l in lines if l.startswith("MODRES")})
    sample = "\n".join(l for l in head
                       if l.startswith("REMARK 210") or l.startswith("REMARK 245")
                       or l.startswith("REMARK 280"))
    ctx = text + "\n" + sample
    rec.update(
        header=hdr, title=title, keywords=kw, expdta=expdta, n_models=nmod,
        n_ssbond=len(ssb), n_link=len(link), het=het, modres=modres,
        is_xray=("X-RAY" in expdta or "ELECTRON" in expdta or "NEUTRON" in expdta),
        is_nmr=("NMR" in expdta),
        fibril=bool(FIBRIL_RE.search(text)),
        lasso=bool(LASSO_RE.search(text)),
        cyclic=bool(CYCLIC_RE.search(text)),
        membrane=bool(MEMB_RE.search(ctx)),
        bound=bool(BOUND_RE.search(title)),
        is_fragment=any(FRAG_RE.search(l) for l in head),
    )
    return rec


def main():
    tg = I.targets()
    recs = []
    for t in tg:
        r = annotate(t["pdb"])
        r.update(n=t["n"], fold=t["fold"], seq=t["seq"],
                 fail18=t["pdb"] in I.FAIL18)
        recs.append(r)
    fail = [r for r in recs if r["fail18"]]
    rest = [r for r in recs if not r["fail18"]]
    flags = ["is_xray", "is_nmr", "fibril", "lasso", "cyclic", "membrane", "bound",
             "is_fragment"]
    tab = {}
    for f in flags:
        a = float(np.mean([r[f] for r in fail])); b = float(np.mean([r[f] for r in rest]))
        tab[f] = dict(fail18=a, other108=b, n_fail=int(sum(r[f] for r in fail)),
                      n_other=int(sum(r[f] for r in rest)))
    for key, fn in (("n_ssbond>0", lambda r: r["n_ssbond"] > 0),
                    ("n_link>0", lambda r: r["n_link"] > 0),
                    ("has_het", lambda r: len(r["het"]) > 0),
                    ("has_modres", lambda r: len(r["modres"]) > 0),
                    ("constrained", lambda r: r["n_ssbond"] > 0 or r["n_link"] > 0
                     or r["lasso"] or r["cyclic"]),
                    ("not_isolated", lambda r: r["fibril"] or r["membrane"] or r["bound"]),
                    ("any_flag", lambda r: (r["fibril"] or r["membrane"] or r["bound"]
                                            or r["lasso"] or r["cyclic"]
                                            or r["n_ssbond"] > 0 or r["n_link"] > 0)),
                    ):
        a = float(np.mean([fn(r) for r in fail])); b = float(np.mean([fn(r) for r in rest]))
        tab[key] = dict(fail18=a, other108=b, n_fail=int(sum(fn(r) for r in fail)),
                        n_other=int(sum(fn(r) for r in rest)))
    out = dict(per_target={r["pdb"]: r for r in recs}, contrast=tab)
    print(json.dumps(tab, indent=1))
    print("\nFAIL18 one-liners:")
    for r in fail:
        print(f"  {r['pdb']} n={r['n']:2d} {r['expdta'][:22]:22s} mdl={r['n_models']:2d} "
              f"ss={r['n_ssbond']} lk={r['n_link']} "
              f"{'FIB' if r['fibril'] else '   '} {'LAS' if r['lasso'] else '   '} "
              f"{'CYC' if r['cyclic'] else '   '} {'MEM' if r['membrane'] else '   '} "
              f"{'BND' if r['bound'] else '   '} | {r['title'][:60]}")
    I.write("fail_headers", out)


if __name__ == "__main__":
    main()
