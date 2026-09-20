import os, sys, json, numpy as np
sys.path.insert(0, r"C:\Users\abena\Protein-Folding-Algorithm")
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"): os.environ.setdefault(v,"1")
from s24 import stats_lib as ST
from s27 import s28_A2_local as A2
sys.path.insert(0, r"C:\Users\abena\Protein-Folding-Algorithm\s30")
import s30_D_gram as G

pdbs=[f[:-4] for f in sorted(os.listdir(G.CACHE)) if f.endswith(".npz")]
rows=[]; names0=None
for pdb in pdbs:
    z=G.load(pdb)
    if names0 is None: names0=z["names"]
    if z["names"]!=names0: continue
    X,keep=G._flat_unit(z["D"]); C0=z["C0"]
    e=A2.remove_rigid(C0-C0.mean(0), C0).ravel(); e=e/np.linalg.norm(e)   # the RADIAL/scale direction
    c=X@e                                                                 # each field's radial cosine
    Gm=X@X.T; tr=float(np.trace(Gm)); w=np.linalg.eigvalsh(Gm)[::-1]
    U,S,Vt=np.linalg.svd(X,full_matrices=False)
    Xr=X-np.outer(c,e)                                                    # radial component removed
    nn=np.linalg.norm(Xr,axis=1); ok=nn>1e-12
    Xr2=Xr[ok]/nn[ok,None]; Gr=Xr2@Xr2.T; wr=np.linalg.eigvalsh(Gr)[::-1]
    uu=z["u"].ravel(); uu=uu/np.linalg.norm(uu)
    rows.append(dict(pdb=pdb, fold=z["fold"], fail18=z["fail18"],
        radial_share=float((c**2).sum()/tr),            # share of the Gram trace along the radial direction
        lam1_share=float(w[0]/tr),
        cos_pc1_radial=float(abs(Vt[0]@e)),             # is the DOMINANT direction the radial one?
        cos_pc2_radial=float(abs(Vt[1]@e)),
        stable=float(tr/w[0]), stable_after=float(wr.sum()/wr[0]),
        lam1_after_share=float(wr[0]/wr.sum()),
        cos_u_radial=float(uu@e),                       # how radial is the direction TO THE NATIVE?
        mean_abs_field_radial=float(np.abs(c).mean())))
folds=ST.pinned_folds([r["pdb"] for r in rows]); fail=np.array([r["fail18"] for r in rows],bool)
def s(k):
    v=np.array([r[k] for r in rows],float)
    c=ST.compare(v,np.zeros(len(v)),folds,seed_parts=("s30Drad",),label=k)
    return v.mean(), np.median(v), c["ci95_fold"], v[fail].mean(), v[~fail].mean()
print("n =",len(rows))
print("%-22s %8s %8s %-20s %8s %8s"%("quantity","mean","median","fold CI","FAIL18","other"))
for k in ["radial_share","lam1_share","cos_pc1_radial","cos_pc2_radial","stable","stable_after",
          "lam1_after_share","cos_u_radial","mean_abs_field_radial"]:
    m,md,ci,f,o=s(k)
    print("%-22s %8.4f %8.4f [%+.3f,%+.3f]   %8.4f %8.4f"%(k,m,md,ci[0],ci[1],f,o))
json.dump(rows, open(sys.argv[1],"w"))
