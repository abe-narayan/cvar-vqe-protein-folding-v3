import os, sys, json, numpy as np
sys.path.insert(0, r"C:\Users\abena\Protein-Folding-Algorithm")
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"): os.environ.setdefault(v,"1")
from s12 import instrument as I
from s24 import stats_lib as ST
rng = np.random.default_rng(3030)
rows=[]
for t in I.targets():
    u=I.load_univ(t["pdb"]); p=I.pool_idx(u); rr=np.asarray(u["rr"][p],float)
    sub=np.asarray(I.shipped_record(t["pdb"])["sub"],int)
    band=np.where(rr<=rr.min()+I.BAND)[0]
    keep=int(np.isin(band,sub).sum())
    # ORACLE best of the score's top-75, and of random 75s drawn from the same 500
    top_best=float(rr[sub].min())
    draws=np.array([rr[rng.choice(len(rr),len(sub),replace=False)].min() for _ in range(200)])
    rows.append(dict(pdb=t["pdb"], fold=int(t["fold"]), fail18=t["pdb"] in I.FAIL18,
                     k=len(rr), nband=int(len(band)), keep=keep,
                     pool_best=float(rr.min()), top_best=top_best,
                     rand_best=float(draws.mean()),
                     delta=float(draws.mean()-top_best)))
json.dump(rows, open(sys.argv[1],"w"))
print("n",len(rows))
