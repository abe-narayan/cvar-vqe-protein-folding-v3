# -*- coding: utf-8 -*-
"""Recompute every headline number in s29/REPORT_S29.md from its artefact.

Standing consequence registered in S29-L48. Prints MATCH/MISMATCH per claim.
"""
import json, glob, io, sys, statistics as st
sys.path.insert(0, '.')
from s24.stats_lib import compare

ok, bad = [], []


def check(label, claimed, actual, tol=5e-4):
    good = abs(claimed - actual) <= tol
    (ok if good else bad).append((label, claimed, actual))
    print('%-58s claim %-10s actual %-10s %s'
          % (label[:58], round(claimed, 4), round(actual, 4), 'MATCH' if good else '*** MISMATCH ***'))


# ---------- lane O ladder ----------
rows = {}
for f in sorted(glob.glob('s29/results/s29_O_chain_rows*.jsonl')):
    for ln in open(f):
        ln = ln.strip()
        if not ln:
            continue
        try:
            r = json.loads(ln)
        except Exception:
            continue
        rows[(r['item'], r['pdb'])] = r
byitem = {}
for (it, pdb), r in rows.items():
    byitem.setdefault(it, {})[pdb] = r


def arm_chain(k):
    d = byitem[k]
    return st.mean(v['rmsd_chain'] for v in d.values()), len(d)


def paired(a, b):
    A, B = byitem[a], byitem[b]
    c = sorted(set(A) & set(B))
    r = compare([A[p]['rmsd_chain'] for p in c], [B[p]['rmsd_chain'] for p in c], names=c, label=a)
    return r


print('=== section 0 / 9 / 11: the ladder (built chain) ===')
check('production, built chain', 3.2105, arm_chain('prod')[0])
check('production, point cloud', 3.0483, st.mean(v['rmsd_cloud'] for v in byitem['prod'].values()))
check('architecture ORACLE ceiling (top-128 prefix avg)', 2.9027, arm_chain('bestm128')[0])
check('best single member, top-128', 2.1435, arm_chain('best1_top128')[0])
check('convex hull, top-128', 1.8538, arm_chain('hull_top128')[0])
check('convex hull, K=500', 1.1235, arm_chain('hull_pool')[0])
r = paired('bestm128', 'prod')
check('ceiling vs production, effect', -0.3079, r['effect'])
check('ceiling vs production, MDE', 0.0982, r['mde'])
r = paired('best1_top128', 'bestm128')
check('best-of-128 vs prefix avg, effect', -0.7592, r['effect'])
check('best-of-128 vs prefix avg, MDE', 0.1815, r['mde'])
check('best-of-128 vs production (the 7-bit contrast)', -1.0670,
      arm_chain('best1_top128')[0] - arm_chain('prod')[0])

# ---------- lane D field survey ----------
print()
print('=== section 9.4: the 21 displacement fields ===')
d = json.load(open('s29/results/s29_D_fields.json'))
check('random-shape reference', 0.1398, d['random_ref_mean_abs'], tol=1e-4)
check('best field cosine (CHAN_DISTPOT)', 0.1128, d['best_mean_cos'])
print('%-58s %s' % ('field count', len(d['fields'])), '(report says 21)')
print('%-58s %s' % ('any field beats the reference', d['any_beats_reference']), '(report says none)')
irs = [v['implied_rmsd_at_best_step'] for v in d['fields'].values() if v.get('implied_rmsd_at_best_step')]
check('best implied RMSD', 3.0289, min(irs))
check('best achievable gain', 0.0195, d['prod_rmsd'] - min(irs))

# ---------- lane B ----------
print()
print('=== section 6: lane B, the set-equality result ===')
brows = []
for f in sorted(glob.glob('s29/results/s29_B_tta_subset_rows.s*of4.jsonl')):
    for ln in open(f):
        ln = ln.strip()
        if ln:
            brows.append(json.loads(ln))
n = len(brows)
print('%-58s %d / %d' % ('non-prefix f-optimal PAIR (report 114/126)',
                         sum(1 for r in brows if not r['pair_is_prefix']), n))
print('%-58s %d / %d' % ('non-prefix f-optimal m=5 SUBSET (report 124/126)',
                         sum(1 for r in brows if not r['greedy_is_prefix']), n))
print('%-58s %d / %d' % ('per-state sort finds optimum (report 12/126)',
                         sum(1 for r in brows if abs(r['perstate_pair_f'] - r['best_pair_f']) < 1e-12), n))
check('mean objective gap, m=5', 0.1499, st.mean(r['greedy_gap_vs_prefix'] for r in brows))
folds = [r['fold'] for r in brows]
names = [r['pdb'] for r in brows]
rr = compare([r['rmsd_greedy'] for r in brows], [r['rmsd_prefix_m'] for r in brows],
             folds=folds, names=names, label='nonprefix')
check('non-prefix choice, effect', 0.0804, rr['effect'])
check('non-prefix choice, MDE', 0.1104, rr['mde'])

# ---------- lane M F1 ----------
print()
print('=== section 4.7: lane M, F1 ===')
try:
    ms = json.load(open('s29/results/s29_M_F1_summary.json'))
    txt = io.open('s29/results/s29_M_F1_fmt.txt', encoding='utf-8').read()
    for needle, lbl in [('+0.0202', 'LOG - SWAPCTL effect +0.0202'),
                        ('+0.0822', 'LOG - PROD effect +0.0822'),
                        ('+0.7295', 'LOGPERM effect +0.7295')]:
        print('%-58s %s' % (lbl, 'FOUND in fmt.txt' if needle in txt else '*** NOT FOUND ***'))
except Exception as e:
    print('F1 artefacts:', e)

# ---------- lane M F2 ----------
print()
print('=== section 12.2: lane M, F2 (the shell profile) ===')
g = json.load(open('s29/results/s29_M_F2_gate.json'))
check('F2 gate, uniform PROD (S12 says 3.078)', 3.0784, g['uniform']['PROD'])
check('F2 gate, uniform ORACLE_PROF (S12 says 2.402)', 2.4023, g['uniform']['ORACLE_PROF'])
check('F2 gate, uniform gap (S12 says 0.676)', 0.6761, g['uniform']['gap'])
check('F2 gate, weighted PROD', 3.0624, g['weighted']['PROD'])
check('F2 gate, weighted ORACLE_PROF', 2.4254, g['weighted']['ORACLE_PROF'])
sup = json.load(open('s29/results/s29_M_F2_supply.json'))['supply']
check('F2 supply, corr(r_hat, r_true)', 0.3129, sup['corr_rhat_rtrue'])
check('F2 supply, corr(r_disto, r_true) [the INCUMBENT]', 0.3660, sup['corr_rdisto_rtrue'])
lam = set(sup['ridge_lambda_per_fold'].values())
print('%-58s %s %s' % ('F2 ridge lambda per fold (report: 1000 on all 5)', sorted(lam),
                       'MATCH' if lam == {1000.0} else '*** MISMATCH ***'))

print()
print('=' * 80)
print('MATCHED: %d    MISMATCHED: %d' % (len(ok), len(bad)))
for lbl, c, a in bad:
    print('  MISMATCH %s: report %s vs artefact %s' % (lbl, c, a))
