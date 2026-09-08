"""s21/rgcheck.py -- THE TWO ATTACKS ON L27 THAT I RAISED AGAINST MY OWN RESULT.

L27 clears a pooled 17-predictor permuted bar.  Two weaknesses in it are mine to test, and I would
rather measure them than wait to be asked.

ATTACK 1 -- THE DIFFICULTY CONTROL CONSUMES THE NATIVE.  L27 partials on `latent_mean`, the mean
RMSD of the whole latent, which reads the native.  I argued that makes the result STRONGER (it
clears a control no deployable pipeline could build).  But that argument cuts both ways: if
difficulty IS partly compactness, an oracle difficulty control may be removing some of the signal --
or, worse, may be leaving structure in that a native-free control would take out.  The honest test
is a NATIVE-FREE difficulty control.

    DECLARED: `pool_spread`, the mean pairwise RMSD within the K=500 retrieval pool.  Purely
    geometric, purely native-free, and a genuine difficulty proxy -- a target whose retrieved
    candidates disagree with each other is a hard target.
    NOT TAKEN: any functional of the shipped score, because that is Workstream D's family and
    partialling my channel on their channel would over-control by construction.

ATTACK 2 -- THE i,i+1 TERM.  The distogram supplies |i-j| >= 2, so my Rg reconstruction adds the
n-1 virtual-bond terms as a CONSTANT 3.8046 A.  The constant is identical for every target, but its
WEIGHT is (n-1)/(n(n-1)/2) = 2/n -- it varies with n, and n is exactly what I residualise on.  If
that interaction is doing the work, the result is a length artefact wearing a compactness hat, which
is the same failure that killed Workstream D's best predictor.

    ABLATION: recompute `rg_disto` with the i,i+1 terms OMITTED entirely.  If rho is unchanged, the
    constant is not driving it.  If rho collapses, L27 is withdrawn.

Both attacks are run at the BINDING control and reported whichever way they land.  Neither is a
robustness garnish: attack 2 has a pre-declared outcome that WITHDRAWS the result.
"""
from __future__ import annotations
import json, os, sys
import numpy as np
for _v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v,"1")
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
if ROOT not in sys.path: sys.path.insert(0,ROOT)
RESULTS=os.path.join(HERE,"results")
from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s21.rgsign import _partial_spearman, _rg2_from_pairs, RG_FAMILY, CA_BOND   # noqa: E402


def run():
    rg={r["pdb"]:r for r in json.load(open(os.path.join(RESULTS,"rgsign.json")))["rows"]}
    lr={r["pdb"]:r for r in json.load(open(os.path.join(RESULTS,"d_lrank.json")))["rows"]}
    P=sorted(set(rg)&set(lr)); rows=[]
    print("targets: %d"%len(P),flush=True)
    for c,pdb in enumerate(P):
        u=I.load_univ(pdb); n=int(u["n"]); p=I.pool_idx(u)
        W=np.asarray(u["W"],float)[p]
        dg=I.distogram(pdb,u["seq"],u["fold"])
        dhat=np.asarray(dg["expected"],float)
        #: ABLATION: no i,i+1 constant at all
        rg_nobond=float(np.sqrt(_rg2_from_pairs(float((dhat**2).sum()),n)))
        #: NATIVE-FREE difficulty: mean pairwise RMSD inside the pool
        sub=W[SD.stable_rng(pdb,"rgcheck").choice(len(W),120,replace=False)]
        Pm=I.pairwise_rmsd(sub)
        pool_spread=float(Pm[np.triu_indices(len(sub),1)].mean())
        r=dict(rg[pdb]); r.update({"rg_nobond":rg_nobond,"pool_spread":pool_spread,
                                   "bond_weight":2.0/n})
        rows.append(r)
        if (c+1)%20==0: print("  %d/%d"%(c+1,len(P)),flush=True)
    ok=len(rows)==len(P) and all(np.isfinite(r["pool_spread"]) for r in rows)
    tmp=os.path.join(RESULTS,"rgcheck.json.tmp")
    json.dump({"rows":rows,"complete":bool(ok),"n_expected":len(P)},open(tmp,"w"))
    os.replace(tmp,os.path.join(RESULTS,"rgcheck.json"))
    report(rows); return rows


def report(rows=None):
    if rows is None:
        rows=json.load(open(os.path.join(RESULTS,"rgcheck.json")))["rows"]
    g=lambda k: np.array([r[k] for r in rows],float)     # noqa: E731
    y=g("ORACLE_rank_pct"); rng=SD.stable_rng("rgcheck","null")
    fam=list(RG_FAMILY)
    print("\nn = %d.  Both attacks at the BINDING control.\n"%len(rows))
    print("  corr(n, bond_weight) = %+.3f   (the term's weight IS a function of length)"
          %np.corrcoef(g("n"),g("bond_weight"))[0,1])
    print("  corr(rg_disto, rg_nobond) = %+.4f\n"%np.corrcoef(g("rg_disto"),g("rg_nobond"))[0,1])
    for lab,ctrls in (("ORACLE difficulty (L27's, latent_mean)",[g("n"),g("latent_mean")]),
                      ("NATIVE-FREE difficulty (pool_spread)", [g("n"),g("pool_spread")]),
                      ("BOTH difficulty controls",             [g("n"),g("latent_mean"),g("pool_spread")])):
        vals={k:_partial_spearman(g(k),y,ctrls) for k in fam}
        vals["rg_nobond [ABLATION]"]=_partial_spearman(g("rg_nobond"),y,ctrls)
        best=max(fam,key=lambda k:abs(vals[k]))
        mx=np.empty(3000)
        for t in range(3000):
            yp=y[rng.permutation(len(y))]
            mx[t]=max(abs(_partial_spearman(g(k),yp,ctrls)) for k in fam)
        bar=float(np.percentile(mx,95))
        print("  CONTROL: %s"%lab)
        for k in fam+["rg_nobond [ABLATION]"]:
            print("    %-24s %+.4f%s"%(k,vals[k]," <- best" if k==best else ""))
        print("    permuted best-of-%d bar %.4f  |  best %.4f  ->  %s\n"
              %(len(fam),bar,abs(vals[best]),"CLEARS" if abs(vals[best])>bar else "NOT DEMONSTRATED"))
    print("PRE-DECLARED: if the ablation collapses, L27 is WITHDRAWN.")


if __name__=="__main__":
    report() if len(sys.argv)>1 and sys.argv[1]=="report" else run()
