# -*- coding: utf-8 -*-
"""Re-collect the 32-cost sweep WITHOUT the fabricated-zero bug.

The first collector used `r.get('pctile') or 0`, which turned a MISSING key into 0.0000.
`native_pctile` is only computed for CA-kind costs (s29_D_cost_audit.py:445-446 run inside the
CA branch), so fifteen chain-kind costs printed a native percentile of exactly 0.0000 and looked
as though the native were each cost's argmin. Fifteen identical values across fifteen different
costs is the same signature that caught lane R's null-input artefact (S30-L19).
"""
import io, json, math, os, sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

COSTS = ['DIS', 'DIS_SURR', 'DIS_MEAN', 'CONTACT_LL', 'DISTPOT', 'CONTACT', 'ENV', 'HP',
         'RG_LAW', 'RG_UNIV', 'EXVOL', 'CAGEO', 'SS_MATCH', 'CONS_POOL', 'DMAP_CONS_POOL',
         'POOLGO_POOL', 'RAMA', 'DSSPHB', 'ELEC', 'TORS_CONS_POOL', 'LEG', 'LEG_steric',
         'LEG_contact', 'LEG_hbond_local', 'LEG_hbond_longrange', 'LEG_coop_helix',
         'LEG_coop_sheet', 'LEG_solvation', 'LEG_electrostatic', 'LEG_aromatic',
         'LEG_torsion', 'LEG_compactness']

NA = None


def dig(d, *k):
    c = d
    for x in k:
        if not isinstance(c, dict) or x not in c:
            return NA
        c = c[x]
    return c


rows, ncomp = [], 0
for c in COSTS:
    fp = 's30/results/s30_D_meter_sw_%s_chain.json' % c
    if not os.path.exists(fp):
        print('MISSING', fp)
        continue
    d = json.load(io.open(fp, encoding='utf-8'))
    rc = dig(d, 'rand_signed_control') or {}
    ct = rc.get('contrast') or {}
    ncomp += dig(d, 'multiplicity', 'comparisons_emitted') or 0
    rows.append(dict(
        cost=c, kind=dig(d, 'kind'),
        ladder=dig(d, 'ladder_rho', 'S28', 'mean'),
        ladder_folds=dig(d, 'ladder_rho', 'S28', 'folds_same_sign'),
        contrast=ct.get('effect'), xmde=ct.get('effect_over_mde'),
        folds=ct.get('folds_same_sign'),
        pctile=dig(d, 'native_pctile', 'mean'),          # None when not computed -- NOT 0
        gate=dig(d, 'verdict', 'gates', 'overall'),
    ))

json.dump(rows, io.open('s30/results/meter_sweep/sweep_chain.json', 'w', encoding='utf-8'),
          indent=1, default=float)

rows.sort(key=lambda r: -(r['xmde'] if r['xmde'] is not None else -9))

print('n costs %d   TOTAL COMPARISONS %d' % (len(rows), ncomp))
# max-over-32 null for |z|: xMDE = z / 2.8016
zmax = max(abs(r['xmde']) for r in rows if r['xmde'] is not None) * 2.8016
print('largest |z| = %.2f;  E[max |z|] over 32 draws under H0 ~ 2.5-2.9 '
      '(i.e. xMDE ~ 0.9-1.04)' % zmax)
print()
print('%-22s %5s %8s %6s %9s %6s %8s %s' %
      ('cost', 'kind', 'ladder', 'lfold', 'contrast', 'xMDE', 'pctile', 'gate'))
for r in rows:
    print('%-22s %5s %8s %6s %9s %6s %8s %s' % (
        r['cost'], r['kind'],
        '%+.4f' % r['ladder'] if r['ladder'] is not None else 'n/a',
        '%s/5' % r['ladder_folds'],
        '%+.4f' % r['contrast'] if r['contrast'] is not None else 'n/a',
        '%+.2f' % r['xmde'] if r['xmde'] is not None else 'n/a',
        '%.4f' % r['pctile'] if r['pctile'] is not None else 'n/a',
        r['gate']))

nca = sum(1 for r in rows if r['pctile'] is not None)
print()
print('native_pctile is computed for %d of %d costs (CA-kind only); the other %d are n/a, '
      'NOT 0.0000' % (nca, len(rows), len(rows) - nca))
