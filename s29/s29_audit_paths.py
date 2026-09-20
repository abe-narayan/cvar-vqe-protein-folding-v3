"""Existence check on every artefact path the S29 ledger cites.

The sprint's own lesson (S29-L41/L42): a path quoted in prose is a claim, not a citation.
This runs the check the contract's "every number carries its artefact path" rule is worth.
"""
import io, os, re, subprocess, collections

s = io.open('s29/LEDGER.md', encoding='utf-8').read()

# Candidate paths: backticked tokens that look like repo paths with a known extension,
# or bare tokens with a directory separator and an extension.
pat = re.compile(r'`([A-Za-z0-9_./\\-]+\.(?:py|json|jsonl|md|csv|pdf|txt|html|npz|pkl))`')
cands = set()
for m in pat.finditer(s):
    p = m.group(1).replace('\\', '/').strip()
    if '/' in p or p.endswith(('.py', '.json', '.jsonl', '.md')):
        cands.add(p)

# Expand shard globs written as name.s*of4.jsonl or *_shard0..N
def exists(p):
    if os.path.exists(p):
        return True
    # tolerate glob-ish citations
    if '*' in p or '..' in p:
        import glob
        base = p.replace('..', '*')
        return bool(glob.glob(base))
    return False

missing, present = [], []
for p in sorted(cands):
    (present if exists(p) else missing).append(p)

print('paths cited in the ledger : %d' % len(cands))
print('present on disk           : %d' % len(present))
print('MISSING                   : %d' % len(missing))
print()

if missing:
    # For each missing path, ask git whether it ever existed -- the S29-L41 check.
    print('%-62s %s' % ('MISSING PATH', 'ever in git history?'))
    print('-' * 90)
    for p in missing:
        try:
            out = subprocess.run(['git', 'log', '--all', '--oneline', '--', p],
                                 capture_output=True, text=True, timeout=25).stdout.strip()
        except Exception as e:
            out = 'check failed: %s' % e
        verdict = ('yes: ' + out.splitlines()[0][:44]) if out else 'NO -- never committed'
        print('%-62s %s' % (p[:62], verdict))

# Which are tracked vs merely on disk?
try:
    tracked = set(subprocess.run(['git', 'ls-files'], capture_output=True, text=True,
                                 timeout=60).stdout.split('\n'))
except Exception:
    tracked = set()
untracked_present = [p for p in present if p not in tracked and '*' not in p and '..' not in p]
print()
print('cited, on disk, but NOT tracked by git: %d' % len(untracked_present))
for p in untracked_present[:25]:
    print('   ', p)
