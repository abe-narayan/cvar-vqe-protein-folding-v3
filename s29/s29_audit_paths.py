# -*- coding: utf-8 -*-
"""Existence check on every artefact path the S29 ledger cites.

The sprint's own lesson (S29-L41/L42): a path quoted in prose is a claim, not a
citation, and it is worth exactly what an existence check on it is worth. This
runs that check over the whole ledger.

A citation resolves if it is (a) a real path on disk, (b) a glob that matches
something, or (c) a basename that matches exactly one file in the tree -- lanes
routinely cite `s29_T_spectra_rows.jsonl` rather than its full path. Anything
left over gets the S29-L41 treatment: `git log --all -- <path>`, to tell "moved"
from "never existed".
"""
import io, os, re, subprocess, glob as _glob

LEDGER = 's29/LEDGER.md'

pat = re.compile(r'`([A-Za-z0-9_./\\-]+\.(?:py|json|jsonl|md|csv|pdf|txt|html|npz|pkl))`')
s = io.open(LEDGER, encoding='utf-8').read()
cands = sorted({m.group(1).replace('\\', '/').strip() for m in pat.finditer(s)})

# index every file in the tree by basename
index = {}
for root, dirs, files in os.walk('.'):
    if '.git' in root:
        continue
    for f in files:
        p = os.path.join(root, f).replace(os.sep, '/')
        index.setdefault(f, []).append(p[2:] if p.startswith('./') else p)


def resolve(p):
    """Return (kind, where) or (None, None)."""
    if os.path.exists(p):
        return 'direct', p
    if '*' in p or '..' in p:
        hits = _glob.glob(p.replace('..', '*'))
        if hits:
            return 'glob', hits[0]
    base = os.path.basename(p)
    hits = index.get(base)
    if hits:
        return 'basename', hits[0]
    # a citation like `name.s*of4.jsonl` written without a directory
    stem = base.split('*')[0]
    loose = [v for k, vs in index.items() if k.startswith(stem) for v in vs]
    if stem and loose:
        return 'prefix', loose[0]
    return None, None


ok, missing = [], []
kinds = {}
for p in cands:
    kind, where = resolve(p)
    if kind:
        ok.append(p)
        kinds[kind] = kinds.get(kind, 0) + 1
    else:
        missing.append(p)

print('paths cited in the ledger : %d' % len(cands))
print('resolve                   : %d  %s' % (len(ok), kinds))
print('DO NOT RESOLVE            : %d' % len(missing))

if missing:
    print()
    print('%-58s %s' % ('UNRESOLVED', 'ever in git history?'))
    print('-' * 90)
    for p in missing:
        try:
            out = subprocess.run(['git', 'log', '--all', '--oneline', '--', p],
                                 capture_output=True, text=True, timeout=25).stdout.strip()
        except Exception as e:
            out = 'check failed: %s' % e
        print('%-58s %s' % (p[:58],
                            ('yes: ' + out.splitlines()[0][:40]) if out else 'NO -- never committed'))

try:
    tracked = set(subprocess.run(['git', 'ls-files'], capture_output=True, text=True,
                                 timeout=60).stdout.split('\n'))
except Exception:
    tracked = set()
untracked = [p for p in ok if os.path.exists(p) and p not in tracked]
print()
print('resolve but are NOT tracked by git: %d' % len(untracked))
for p in untracked[:20]:
    print('   ', p)
