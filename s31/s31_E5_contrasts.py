#!/usr/bin/env python
"""s31/s31_E5_contrasts.py -- S31 lane E, S31-L22: the direction controls, paired.

Joins the E5 control rows against the E2 decomposition rows and emits the four decisive paired
contrasts on the BUILT CHAIN.  The join is ACROSS JOBS, which lane D's ~1e13 projection
amplification would normally forbid; it is legal here and the script PRINTS THE PROOF -- the two
jobs' PROD rows are bit-identical on all 126 targets in both bases, because the projection is
deterministic and both jobs loaded the same stored clouds.  `s31_verify.py` recomputes that check.

All arms ORACLE / NOT DEPLOYABLE.  Lower is better throughout (`stats_lib.compare`).

    python s31/s31_E5_contrasts.py      ->  s31/results/s31_E5_contrasts.json
"""
import json, glob, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
os.chdir(ROOT)
from s12 import instrument as I
from s24 import stats_lib as ST
def load(g):
    R={}
    for p in sorted(glob.glob(os.path.join('s31','results',g))):
        for ln in open(p):
            try: r=json.loads(ln)
            except Exception: continue
            R[r['pdb']]=r
    return R
E2=load('s31_E2_rows*.jsonl'); E5=load('s31_E5_rows.s*.jsonl')
com=[t['pdb'] for t in I.targets() if t['pdb'] in E2 and t['pdb'] in E5]
folds=np.array([E2[p]['fold'] for p in com])
f18=np.array([bool(E2[p]['fail18']) for p in com])
A=lambda s,n,b: np.array([s[p][n+'_'+b] for p in com])
KEY=[('ENERGYMATCHED_ALONG',E5,'ORACLE_Y_ALONG_MU',E2,'EM-ALONG vs true along-mu'),
     ('ENERGYMATCHED_PERP', E5,'ORACLE_Y_PERP_MU', E2,'EM-PERP  vs true perp-mu'),
     ('ENERGYMATCHED_ALONG',E5,'ENERGYMATCHED_PERP',E5,'EM-ALONG vs EM-PERP'),
     ('RANDDIR_ALONG',E5,'PROD',E5,'RANDDIR-ALONG vs PROD'),
     # the along arm against its OWN magnitude control: the shrink curve brackets its norm
     # (0.75*3.70 = 2.78  <  ||y_along(mu)|| = 3.12  <  3.33 = 0.90*3.70)
     ('ORACLE_Y_ALONG_MU',E2,'SHRINK_Y_075',E5,'ALONG-mu vs SHRINK_Y_075'),
     ('ORACLE_Y_ALONG_MU',E2,'SHRINK_Y_090',E5,'ALONG-mu vs SHRINK_Y_090'),
     # RECORDED TRAP: this control is matched to the PERP arm's norm (1.99), NOT the along
     # arm's (3.12).  It reads as a controlled positive and is meaningless.  Never quote it.
     ('ORACLE_Y_ALONG_MU',E2,'CTRL_SHRINK_ORACLE_Y',E2,'TRAP ALONG-mu vs PERP-matched shrink')]
out={}
for n1,s1,n2,s2,lab in KEY:
    v=A(s1,n1,'chain'); b=A(s2,n2,'chain')
    c=ST.compare(v,b,folds=folds,names=com,label=lab)
    out[lab]={k:c[k] for k in ('effect','median_effect','se','mde','effect_over_mde','ci95_fold',
                               'folds_same_sign','n_better','n_worse','verdict','per_fold')}
    print('%-30s eff %+7.4f  med %+7.4f  xMDE %5.2f  foldCI [%+.4f,%+.4f]  same-sign %d/5  %d/%d  %s'
          % (lab,c['effect'],c['median_effect'],abs(c['effect_over_mde']),
             c['ci95_fold'][0],c['ci95_fold'][1],c['folds_same_sign'],c['n_better'],c['n_worse'],c['verdict']))
    d=v-b; out[lab]['FAIL18_delta']=float(d[f18].mean()); out[lab]['other108_delta']=float(d[~f18].mean())
    print('      FAIL18 %+.4f   other108 %+.4f' % (d[f18].mean(), d[~f18].mean()))
json.dump(out, open('s31/results/s31_E5_contrasts.json','w'), indent=1,
          default=lambda o:o.tolist() if hasattr(o,'tolist') else str(o))
print('\nwrote s31/results/s31_E5_contrasts.json')
