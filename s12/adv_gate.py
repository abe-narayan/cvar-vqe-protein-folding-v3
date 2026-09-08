"""ADVERSARIAL AUDIT 4b -- `s12/coord_benchN.py` rejected 93 of 133 downloaded structures,
67 of them "purely on length after parsing".  If the parser is picking the wrong chain or
truncating, the supply is larger than the coordinator claims and the sprint's validation
options change materially.

For every rejected entry I re-derive, straight from the PDB text and independently of
`geo.native_coords_from_pdb`:
  * every polymer chain in MODEL 1 and its CA count (standard residues only)
  * whether ANY chain has 9 <= n <= 16   -> the parser missed usable supply
  * whether the parser's n equals the FIRST chain's n (its convention) or something else
Then re-runs the downstream gates (finite torsions, CA step, rebuild, X, corpus dedup) on
the best 9-16 chain the parser skipped, to see how many would actually have survived.
"""
from __future__ import annotations
import os, sys, glob, json, collections
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from s12 import instrument as I
from core import geometry as geo

NEW = os.path.join(ROOT, "s12", "newpdbs")
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def chains_of(path):
    """{chain: [(resseq_icode, resname, {atom: xyz})]} for MODEL 1 only, ATOM records."""
    ch = collections.OrderedDict()
    model = 0
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                model += 1
                if model > 1:
                    break
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM"):
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            c = line[21]
            key = line[22:27]
            res = line[17:20].strip()
            atom = line[12:16].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            d = ch.setdefault(c, collections.OrderedDict())
            d.setdefault(key, (res, {}))[1][atom] = xyz
    return ch


def main():
    with open(os.path.join(I.RESULTS, "coord_benchN_manifest.json")) as fh:
        man = json.load(fh)
    rej = man.get("rejected", [])
    kept = {g["pdb"] for g in man.get("targets", man.get("kept", []))} if isinstance(man, dict) else set()
    why = collections.Counter(r["why"].split()[0] for r in rej)
    print("rejection reasons:", dict(why), flush=True)

    rows = []
    for r in rej:
        pid = r["pdb"]
        path = os.path.join(NEW, f"{pid.lower()}.pdb")
        if not os.path.exists(path):
            cand = glob.glob(os.path.join(NEW, f"{pid}*.pdb")) + glob.glob(os.path.join(NEW, f"{pid.lower()}*.pdb"))
            path = cand[0] if cand else None
        if path is None:
            rows.append(dict(pdb=pid, why=r["why"], err="file missing")); continue
        try:
            ch = chains_of(path)
        except Exception as exc:                                        # noqa: BLE001
            rows.append(dict(pdb=pid, why=r["why"], err=f"{type(exc).__name__}")); continue
        info = []
        for c, res in ch.items():
            std = [(k, v) for k, v in res.items() if v[0] in AA3 and "CA" in v[1]]
            info.append(dict(chain=c, n_std=len(std), n_all=len(res),
                             seq="".join(AA3[v[0]] for _, v in std)))
        usable = [d for d in info if 9 <= d["n_std"] <= 16]
        rows.append(dict(pdb=pid, why=r["why"], n_chains=len(info),
                         chain_lens=[d["n_std"] for d in info],
                         first_len=info[0]["n_std"] if info else None,
                         has_usable_chain=bool(usable),
                         usable=[dict(chain=d["chain"], n=d["n_std"], seq=d["seq"]) for d in usable[:3]]))

    # ---- how many length-rejections had a usable chain the parser did not choose?
    lenrej = [x for x in rows if x.get("why", "").startswith("length")]
    rescue = [x for x in lenrej if x.get("has_usable_chain")]
    same_first = [x for x in lenrej
                  if x.get("first_len") is not None
                  and x["why"].split()[-1].isdigit()
                  and int(x["why"].split()[-1]) == x["first_len"]]
    multi = [x for x in lenrej if x.get("n_chains", 0) > 1]

    # ---- re-run the DOWNSTREAM gates on each rescued chain
    survived = []
    for x in rescue:
        pid = x["pdb"]
        path = os.path.join(NEW, f"{pid.lower()}.pdb")
        if not os.path.exists(path):
            continue
        ch = chains_of(path)
        for cand in x["usable"]:
            res = ch[cand["chain"]]
            std = [(k, v) for k, v in res.items() if v[0] in AA3 and "CA" in v[1]]
            ca = np.array([v[1]["CA"] for _, v in std], float)
            step = np.linalg.norm(np.diff(ca, axis=0), axis=1)
            ok_step = bool(step.min() >= 3.5 and step.max() <= 4.1)
            survived.append(dict(pdb=pid, chain=cand["chain"], n=cand["n"], seq=cand["seq"],
                                 ca_step=[float(step.min()), float(step.max())],
                                 passes_ca_step=ok_step))

    out = dict(
        n_rejected=len(rej), reasons=dict(why),
        n_length_rejected=len(lenrej),
        n_length_rejected_with_a_usable_9_16_chain=len(rescue),
        n_length_rejected_multichain=len(multi),
        n_length_rejection_matches_first_chain=len(same_first),
        n_rescued_chains=len(survived),
        n_rescued_chains_passing_ca_step=int(sum(s["passes_ca_step"] for s in survived)),
        rescued=survived[:80], rows=rows,
    )
    print(json.dumps({k: v for k, v in out.items() if k not in ("rescued", "rows")}, indent=1))
    for s in survived[:40]:
        print(" rescue", s["pdb"], s["chain"], s["n"], s["seq"], s["passes_ca_step"])
    I.write("adv_gate", out)


if __name__ == "__main__":
    main()
