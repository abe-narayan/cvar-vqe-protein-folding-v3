/* s25/resultslab/site/test_site.js -- the browser's geometry against Python's.
 *
 *   node s25/resultslab/site/test_site.js  [results/site]
 *
 * The page displays the RMSD Python computed from the exported bytes. This test asserts the
 * INDEPENDENT browser-side recomputation (Horn quaternion superposition, proper rotations
 * only) agrees with it, so the tripwire in `rmsdBar` can never be the thing that is wrong.
 */
'use strict';
const path = require('path');
const fs = require('fs');

const siteDir = process.argv[2] || path.join(__dirname, '..', '_selftest', 'site');
global.window = {addEventListener() {}};
require(path.resolve(siteDir, 'data.js'));
require(path.resolve(siteDir, 'structures.js'));
const APP = require(path.resolve(siteDir, 'app.js'));

const D = window.RESULTS_DATA, S = window.STRUCTURES;
let n = 0, worst = 0, worstCell = null;
for (const r of D.records) {
  if (r.status !== 'present') continue;
  const nat = APP.parsePDB(S[`${r.target_id}__native__native_ca`]).ca;
  if (!nat.length) { console.error('MISSING native for', r.target_id); process.exit(1); }
  // BOTH bases: the key carries the basis, so a mix-up is impossible rather than unlikely
  for (const [basis, value] of [[r.basis, r.rmsd], [r.basis_secondary, r.rmsd_secondary]]) {
    const txt = S[`${r.target_id}__${r.configuration}__${basis}`];
    if (!txt) { console.error('MISSING structure', r.target_id, r.configuration, basis); process.exit(1); }
    const parsed = APP.parsePDB(txt);
    if (parsed.ca.length !== r.length) {
      console.error('LENGTH mismatch', r.target_id, r.configuration, basis, parsed.ca.length, r.length);
      process.exit(1);
    }
    if (parsed.meta.basis !== basis) {
      console.error('BASIS mismatch: key says', basis, 'header says', parsed.meta.basis);
      process.exit(1);
    }
    if (!parsed.meta.provenance || !parsed.meta.git_commit || !parsed.meta.module_sha) {
      console.error('PROVENANCE missing from header of', r.target_id, r.configuration, basis);
      process.exit(1);
    }
    const {rmsd} = APP.superpose(parsed.ca, nat);
    const d = Math.abs(rmsd - value);
    if (d > worst) { worst = d; worstCell = `${r.target_id}/${r.configuration}/${basis}`; }
    n++;
  }
}
const TOL = 1e-6;
console.log(`  browser vs python RMSD over ${n} cells: worst |Δ| = ${worst.toExponential(3)} Å  (${worstCell})`);
if (worst > TOL) { console.error(`  FAIL: exceeds ${TOL}`); process.exit(1); }

/* a reflection must NOT be used: mirroring a structure must not lower its RMSD to zero */
const t0 = D.records.find(r => r.status === 'present');
const nat = APP.parsePDB(S[`${t0.target_id}__native__native_ca`]).ca;
const mirrored = nat.map(p => [-p[0], p[1], p[2]]);
const {rmsd: mr} = APP.superpose(mirrored, nat);
if (!(mr > 0.5)) { console.error('  FAIL: a reflected copy superposed onto the native (RMSD ' + mr + ')'); process.exit(1); }
console.log(`  reflection check: mirrored native superposes at ${mr.toFixed(4)} Å, not 0, so proper rotations only`);

/* every leaderboard row the page will render must carry the keys it reads */
const need = ['configuration', 'n_present', 'n_targets', 'mean', 'median', 'sd', 'verdict',
              'mean_secondary', 'basis_secondary', 'pool_gate', 'difficulty_gate', 'corr_with_pool_best'];
for (const row of D.leaderboard) for (const k of need) {
  if (!(k in row)) { console.error('  FAIL: leaderboard row missing', k); process.exit(1); }
}
console.log('  site tests PASS');
