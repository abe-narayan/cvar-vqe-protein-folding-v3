# -*- coding: utf-8 -*-
"""Recompute every headline number of Sprint 30 from its own artefact.

Standing consequence carried from S29: a verdict string is a claim about a
computation and is worth exactly what reading the computation is worth. S30 found
that the hard way twice -- once when an S29 verdict tested the wrong statistic
against the wrong null (S30-L5), and once when an entry's headline was conditioned
on its own numerator (S30-L23).

Usage:  python s30/s30_verify.py
"""
import io, json, os, sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

OK, BAD, MISSING = [], [], []


def load(path):
    if not os.path.exists(path):
        MISSING.append(path)
        return None
    try:
        return json.load(io.open(path, encoding='utf-8'))
    except Exception as e:
        MISSING.append('%s (%s)' % (path, e))
        return None


def dig(obj, *keys, **kw):
    """Walk nested dicts; return default if any key is absent."""
    cur = obj
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return kw.get('default')
        cur = cur[k]
    return cur


def check(label, claimed, actual, tol=5e-4):
    if actual is None:
        MISSING.append(label)
        print('%-56s claim %-10s  ARTEFACT KEY NOT FOUND' % (label[:56], claimed))
        return
    good = abs(float(claimed) - float(actual)) <= tol
    (OK if good else BAD).append((label, claimed, actual))
    print('%-56s claim %-10s actual %-10s %s'
          % (label[:56], round(float(claimed), 4), round(float(actual), 4),
             'MATCH' if good else '*** MISMATCH ***'))


def show(label, value, note=''):
    print('%-56s %-22s %s' % (label[:56], str(value)[:22], note))


print('=' * 92)
print('SPRINT 30 -- recomputing headline numbers from artefacts')
print('=' * 92)

# ---------------------------------------------------------------- lane D: the Gram
print()
print('--- lane D, the field Gram (S30-L21): the combination question ---')
g = load('s30/results/s30_D_gram.json')
if g:
    check('Gram stable rank', 2.0572, dig(g, 'gram_mean', 'stable_rank'), tol=5e-3)
    check('per-target stable rank', 1.6811, dig(g, 'per_target_rank', 'stable'), tol=5e-3)
    check('ORACLE GLOBAL combination rho', 0.1693, dig(g, 'global_weights', 'oracle_global_rho'))
    check('leave-fold-out combination rho', 0.0124, dig(g, 'global_weights', 'lfo_rho'))
    check('best SINGLE field rho (the floor)', 0.1214,
          dig(g, 'global_weights', 'lfo_best_single_rho'))
    check('ORACLE PER-TARGET rho (the trap)', 0.9491, dig(g, 'oracle_per_target', 'rho'))
    check('  its matched random 21-dim control', 0.8095,
          dig(g, 'oracle_per_target', 'ctrl_full_dim_rho'))
    show('  -> LFO is BELOW the best single field',
         dig(g, 'global_weights', 'lfo_rho', default=9) <
         dig(g, 'global_weights', 'lfo_best_single_rho', default=0),
         'the 21 fields contain no combinable information')

# ---------------------------------------------------------------- lane D: the radial test
print()
print("--- lane D, lane L's radial prediction (S30-L22): per-target rows ---")
r = load('s30/results/s30_D_radial.json')
if isinstance(r, list) and r:
    import statistics as st
    f18 = [w for w in r if w.get('fail18')]
    # lane D reported both: mean 0.947, median 0.984 -- compare like with like
    check('cos(dominant principal direction, radial), MEAN', 0.9472,
          st.mean([w['cos_pc1_radial'] for w in r]), tol=2e-2)
    check('cos(dominant principal direction, radial), MEDIAN', 0.9840,
          st.median([w['cos_pc1_radial'] for w in r]), tol=2e-2)
    check('cos(direction to native, radial), all, mean', -0.0675,
          st.mean([w['cos_u_radial'] for w in r]), tol=2e-2)
    check('cos(direction to native, radial), FAIL18, mean', -0.2524,
          st.mean([w['cos_u_radial'] for w in f18]) if f18 else None, tol=2e-2)
    show('  n rows / n FAIL18', '%d / %d' % (len(r), len(f18)),
         'the library spends its rank on a direction anti-aligned with the answer')

# ---------------------------------------------------------------- lane T: bits
print()
print('--- lane T, the ORACLE ladder behind the value-of-a-bit law (S30-L14) ---')
t = load('s30/results/s30_T_bits.json')
if isinstance(t, dict) and isinstance(t.get('rows'), list):
    rows = t['rows']
    show('  ladder rows', len(rows), 'fit reported at R^2 0.9983, a = 1.3312, gamma = 3.1636')
    if rows and isinstance(rows[0], dict):
        show('  row keys', sorted(rows[0].keys())[:8], '')

# ---------------------------------------------------------------- lane X
print()
print('--- lane X, the set-mean decomposition (S30-L20) ---')
x = load('s30/results/s30_X_typicalgood.json')
if isinstance(x, dict):
    pa = x.get('per_arm')
    if isinstance(pa, dict):
        for arm in ['T0_helix', 'T1_blind', 'T2_restype', 'T3_pool']:
            a = pa.get(arm)
            if isinstance(a, dict):
                B = a.get('B_endpoint', a.get('B', a.get('endpoint')))
                show('  %s' % arm,
                     'set_mean %s  B %s  S %s' % (round(a.get('set_mean', 0), 4),
                                                  'MISSING' if B is None else round(B, 4),
                                                  round(a.get('S', 0), 4)), '')
    show('  MOST_CONCENTRATED_SOURCE', str(x.get('MOST_CONCENTRATED_SOURCE'))[:60], '')
    show('  SECONDARY corr(S,B) within target', str(x.get('SECONDARY_within_target_corr_S_B'))[:60],
         'lane X: the strong negative is a T0_helix artefact')

# ---------------------------------------------------------------- lane R
print()
print('--- lane R, the recognition verdict (S30-L19) ---')
rv = load('s30/results/s30_R_verdict.json')
if isinstance(rv, dict):
    show('  n targets', rv.get('n'), '')
    d1 = rv.get('D1_locality')
    if isinstance(d1, dict):
        show('  D1 locality', str(d1)[:110], 'local dR2 -0.089 vs global +0.600')
    show('  D3 anchor confound', str(rv.get('D3_anchor_confound'))[:80],
         'threshold 0.5; lane R reports median +0.259, did not bind')

# ---------------------------------------------------------------- lane F
print()
print('--- lane F, the stage-gap decomposition (S30-L2) ---')
f = load('s30/results/s30_F_stagegap.json')
if isinstance(f, dict):
    show('  n / n_fail18', '%s / %s' % (f.get('n'), f.get('n_fail18')), '')
    gaps = f.get('gaps')
    if isinstance(gaps, dict):
        show('  gaps keys', list(gaps.keys())[:8], '')
    show('  F1a (pool-limited?)', str(f.get('F1a'))[:90],
         'ORACLE best of FAIL18 pools 2.2842 -- NOT pool-limited')
    show('  F1c', str(f.get('F1c'))[:90],
         'WITHDRAWN as an effect estimate by S30-L23 -- the stratum is the outcome')

# ------------------------------------------------- lane D: the meter, and the charter's anchors
print()
print("--- lane D, the cost/RMSD meter (S30-L24): the charter's four anchors ---")
mc = load('s30/results/s30_D_meter_DIS_chain.json')
ma = load('s30/results/s30_D_meter_DIS_ca.json')
if mc:
    check("charter anchor: ladder rho (chain, S28)", -0.40,
          dig(mc, 'ladder_rho', 'S28', 'mean'), tol=1e-2)
    check("charter anchor: ORACLE preference (chain)", 0.07,
          dig(mc, 'rand_signed_control', 'pref_circ_best'), tol=5e-3)
    check("production RMSD, built chain", 3.2071, dig(mc, 'mean_rmsd', 'PROD'), tol=5e-3)
if ma:
    check("charter anchor: cosine (ca)", -0.03, dig(ma, 'cosine', 'mean'), tol=5e-3)
    check("charter anchor: native percentile", 0.369,
          dig(ma, 'native_pctile', 'mean'), tol=2e-3)
    check("production RMSD, CA point cloud", 3.0483, dig(ma, 'mean_rmsd', 'PROD'), tol=5e-3)

# The S30-L24 withdrawal: 8 draws, and the two bases split.
print()
print('--- S30-L24: the 8-draw random-signed control, and why the bases split ---')
for nm, o, claimed_eff, claimed_x in (('chain', mc, 0.0357, 0.5636), ('ca', ma, 0.1716, 1.8415)):
    if not o:
        continue
    rc = dig(o, 'rand_signed_control') or {}
    ct = rc.get('contrast') or {}
    check('  %-5s preference contrast (8 draws)' % nm, claimed_eff, ct.get('effect'), tol=1e-3)
    check('  %-5s   as a multiple of its own MDE' % nm, claimed_x, ct.get('effect_over_mde'), tol=2e-3)
    show('  %-5s   n_draws / draw sd' % nm,
         '%s / %s' % (rc.get('n_draws'), round(rc.get('pref_rand_signed_draw_sd') or 0, 4)),
         'single-draw range %s' % (rc.get('single_draw_range'),))
    sd_rel = ((rc.get('pref_rand_signed_draw_sd') or 0)
              / (rc.get('pref_rand_signed_mean') or 1))
    show('  %-5s   RELATIVE draw noise' % nm, '%.1f%%' % (100 * sd_rel),
         'the chain is the noisy basis -- the control is a rare event there')

# ---- the sign-convention trap, asserted so it stays known rather than lurking ----
print()
print('--- the sign-convention trap in this artefact (KNOWN, asserted, not a defect) ---')
if ma:
    rc = dig(ma, 'rand_signed_control') or {}
    ct = rc.get('contrast') or {}
    gate = dig(ma, 'verdict', 'gates', 'rand_signed')
    vstr = str(ct.get('verdict'))
    # s24.stats_lib.compare is LOWER-IS-BETTER (d = a - b, negative = a better) because its
    # native statistic is RMSD.  A PREFERENCE RATE is HIGHER-IS-BETTER, so for this one
    # statistic its verdict string is inverted.  Lane D's gate handles it correctly
    # (s30_D_meter.py:451, `PASS if effect > 0`).  Anyone quoting `.verdict` instead of the
    # gate reports the OPPOSITE of the truth.  No S30 document does; this keeps it that way.
    ok = (vstr.startswith('WORSE') and gate == 'PASS' and (ct.get('effect') or 0) > 0)
    (OK if ok else BAD).append(('sign-convention trap present and gated', 'WORSE/PASS',
                                '%s / %s' % (vstr[:12], gate)))
    print('%-56s %s' % ('  contrast.verdict says', vstr[:40]))
    print('%-56s %s' % ('  verdict.gates.rand_signed says', gate))
    print('%-56s %s' % ('  -> QUOTE THE GATE, NOT .verdict',
                        'CONSISTENT with the known trap' if ok else '*** trap changed shape ***'))

# ------------------------------------------- the ENDPOINT itself, and its reproducibility (D10)
print()
print('--- the endpoint: production, and how well it reproduces ---')
import statistics as _st


def _rows(path, key, **flt):
    if not os.path.exists(path):
        MISSING.append(path)
        return None
    out = []
    for line in io.open(path, encoding='utf-8'):
        line = line.strip()
        if not line:
            continue
        try:
            r = json.loads(line)
        except Exception:
            continue
        if all(r.get(k) == val for k, val in flt.items()):
            out.append(r[key])
    return out or None


s29p = _rows('s29/results/s29_O_chain_rows.jsonl', 'rmsd_chain', item='prod')
s29c = _rows('s29/results/s29_O_chain_rows.jsonl', 'rmsd_cloud', item='prod')
s27p = _rows('s27/results/chain_rows.jsonl', 'rmsd_chain', config='DIS')
s27c = _rows('s27/results/chain_rows.jsonl', 'rmsd_cloud', config='DIS')

if s29p:
    check('production built chain (S29 prod row, CANONICAL)', 3.2105, _st.mean(s29p))
    check('production CA cloud   (S29 prod row)', 3.0483, _st.mean(s29c))
    show('  n targets', len(s29p), '')
if s27p:
    check('production built chain (S27 DIS config)', 3.2126, _st.mean(s27p))
    check('production CA cloud   (S27 DIS config)', 3.0483, _st.mean(s27c))
if s29p and s27p:
    dchain = abs(_st.mean(s29p) - _st.mean(s27p))
    dcloud = abs(_st.mean(s29c) - _st.mean(s27c))
    show('  chain disagreement between records', '%.4f A' % dchain,
         'the projection is multi-start and NOT seed-pinned')
    show('  cloud disagreement between records', '%.4f A' % dcloud,
         'the cloud is EXACT; only the projection is stochastic')
    ok = dcloud < 1e-4 < dchain
    (OK if ok else BAD).append(('cloud exact / chain stochastic', 'cloud==, chain!=',
                                '%.4f / %.4f' % (dcloud, dchain)))
    show('  -> any claim below ~0.01 A on the chain', 'IS INSIDE THE NOISE',
         'E2 at 0.0221 is only ~2x the full 0.0107 spread')

# ------------------------------------------- every path the ledger claims to have written (S30-L0)
print()
print('--- every path the ledger claims to have written (the S30-L0 failure) ---')
claimed = ['s30/BRIEF.md', 's30/LEDGER.md', 's30/STATE.md', 's30/THEORY.md',
           's30/S30_CONTRACT.md', 's30/REPORT_S30.md', 's30/s30_D_meter.py',
           's30/AUDIT_V.md', 's30/QUANTUM_W.md']
for c in claimed:
    good = os.path.exists(c)
    (OK if good else BAD).append(('path exists: %s' % c, 'exists', 'yes' if good else 'MISSING'))
    print('%-56s %s' % ('  ' + c, 'exists' if good else '*** MISSING ***'))

print()
print('=' * 92)
print('MATCHED: %d    MISMATCHED: %d    KEYS/FILES NOT FOUND: %d' % (len(OK), len(BAD), len(MISSING)))
for lbl, c, a in BAD:
    print('  MISMATCH  %s: recorded %s vs artefact %s' % (lbl, c, a))
for m in MISSING[:12]:
    print('  NOT FOUND %s' % m)
print('=' * 92)
