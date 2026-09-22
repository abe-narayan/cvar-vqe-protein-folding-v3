"""L2: what experimental method determined each of the 126 targets?

Decides the new-observable question. For an NMR-solved peptide the deposited restraints ARE
what fixed the coordinates, so using them is the native through a different door, not a new
channel.  For an X-ray structure the same is true of the density.
"""
import sys, json, collections
import urllib.request

import os as _os; sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))))
import core.pipeline as P

ids = [t.pdb.upper() for t in P.manifest("tuning126")]
print(f"{len(ids)} targets, {len(set(ids))} unique")

Q = """query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    exptl { method }
    rcsb_entry_info { resolution_combined deposited_polymer_monomer_count }
    pdbx_nmr_ensemble { conformers_calculated_total_number }
    rcsb_entry_container_identifiers { related_emdb_ids }
  }
}"""

out = {}
B = 40
for i in range(0, len(ids), B):
    chunk = ids[i:i + B]
    body = json.dumps({"query": Q, "variables": {"ids": chunk}}).encode()
    req = urllib.request.Request(
        "https://data.rcsb.org/graphql", data=body,
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        d = json.load(r)
    for e in d["data"]["entries"] or []:
        if e is None:
            continue
        meth = [m["method"] for m in (e.get("exptl") or [])]
        out[e["rcsb_id"]] = dict(
            method=";".join(meth) or "?",
            res=(e.get("rcsb_entry_info") or {}).get("resolution_combined"),
            nres=(e.get("rcsb_entry_info") or {}).get("deposited_polymer_monomer_count"),
            nconf=(e.get("pdbx_nmr_ensemble") or {}).get(
                "conformers_calculated_total_number"),
        )

miss = [i for i in ids if i not in out]
cnt = collections.Counter(v["method"] for v in out.values())
print("\nEXPERIMENTAL METHOD over tuning126:")
for k, v in cnt.most_common():
    print(f"  {v:4d}  ({100*v/len(ids):5.1f}%)  {k}")
if miss:
    print(f"  {len(miss):4d}  NOT RETURNED BY RCSB: {miss[:10]}")

xr = [k for k, v in out.items() if "X-RAY" in v["method"].upper()]
nmr = [k for k, v in out.items() if "NMR" in v["method"].upper()]
print(f"\nNMR-determined: {len(nmr)}/{len(ids)} = {100*len(nmr)/len(ids):.1f}%")
print(f"X-ray:          {len(xr)}/{len(ids)} = {100*len(xr)/len(ids):.1f}%")

json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "expmethod.json", "w"), indent=1)
print("\nwrote", sys.argv[1] if len(sys.argv) > 1 else "expmethod.json")
