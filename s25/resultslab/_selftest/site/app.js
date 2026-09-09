/* s25 results laboratory — browse flow + 3D protein viewer.
 *
 * Data contract (written by s25/resultslab/site.py):
 *   window.RESULTS_DATA = {meta, targets[], leaderboard[], records[]}
 *   window.STRUCTURES   = {"<T-number>__<configuration>__<basis>": "<pdb text>", ...}
 *
 * THE BASIS IS PART OF EVERY STRUCTURE KEY. That is not tidiness: s25 L4 happened because a
 * number travelled without its basis. Here the viewer physically cannot load a structure on
 * one basis while printing the RMSD of the other, and every RMSD on screen is rendered with
 * its basis attached.
 *
 * WHAT THIS PAGE IS FOR. The headline of this project is a NEGATIVE result with a
 * well-characterised mechanism. The page is not a victory lap; it is built so a reader can
 * CHECK it -- every number traces to a file, every file to a commit, every basis is stated.
 *
 * DESIGN NOTES THAT MATTER SCIENTIFICALLY
 *  - The RMSD shown is ALWAYS the one Python computed from the bytes of the exported file
 *    (record.rmsd, cross-checked against the file's own REMARK 999 RMSD_CA header). The
 *    browser recomputes it independently as a tripwire and flags any disagreement > 0.01 A
 *    rather than quietly preferring one of them.
 *  - The BUILT CHAIN is the production result. The POINT CLOUD is shown beside it, always,
 *    labelled a non-physical intermediate -- its virtual Ca-Ca bonds are 22% contracted.
 *  - Structure files hold predictions in their OWN frame. Superposition onto the native is
 *    done here, at draw time, by Horn's quaternion method (proper rotations only — a
 *    reflection would make an L chain look like a D chain and lower the RMSD illegitimately).
 *  - Camera and representation are identical for every configuration of a target: both
 *    viewers frame the NATIVE, so structures are visually comparable without relearning.
 */
'use strict';

const D = window.RESULTS_DATA || {meta: {}, targets: [], leaderboard: [], records: []};
const S = window.STRUCTURES || {};
const HAS3D = typeof window.$3Dmol !== 'undefined';

const byTarget = new Map();
const byConfig = new Map();
D.records.forEach(r => {
  if (!byTarget.has(r.target_id)) byTarget.set(r.target_id, new Map());
  byTarget.get(r.target_id).set(r.configuration, r);
  if (!byConfig.has(r.configuration)) byConfig.set(r.configuration, []);
  byConfig.get(r.configuration).push(r);
});
const CONFIGS = D.meta.configurations || [...byConfig.keys()];
const BASELINE = D.meta.baseline;
const B1 = D.meta.basis;              // built_chain_bb  -- the production result
const B2 = D.meta.basis_secondary;    // point_cloud_ca  -- non-physical intermediate
const BNAT = 'native_ca';
const skey = (tid, cfg, basis) => `${tid}__${cfg}__${basis}`;
const basisLabel = b => b === B1 ? 'built chain (production)'
                      : b === B2 ? 'point cloud (non-physical intermediate)'
                      : b;

/* ------------------------------------------------------------------ small helpers */
const el = (t, a = {}, kids = []) => {
  const n = document.createElement(t);
  for (const [k, v] of Object.entries(a)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'text') n.textContent = v;
    else if (k === 'html') n.innerHTML = v;
    else if (k === 'cls') n.className = v;
    else n.setAttribute(k, v === true ? '' : v);
  }
  (Array.isArray(kids) ? kids : [kids]).forEach(c => c && n.append(c));
  return n;
};
const num = (x, d = 4) => (x === null || x === undefined || Number.isNaN(x)) ? '—' : Number(x).toFixed(d);
const signed = (x, d = 4) => (x === null || x === undefined) ? '—' : (x >= 0 ? '+' : '') + Number(x).toFixed(d);
// A real configuration name is a "+"-joined combination of Hamiltonians; a TEST FIXTURE
// label starts with "_" and is shown verbatim so it can never read as a method.
const pretty = c => c.startsWith('_') ? c : c.replace(/_/g, ' + ');

/* ------------------------------------------------------------------ geometry (JS) */
function parsePDB(text) {
  const meta = {}, ca = [], names = [], all = [];
  for (const line of text.split('\n')) {
    const m = /^REMARK 999 ([A-Z0-9_]+)\s(.*)$/.exec(line);
    if (m) { meta[m[1].toLowerCase()] = m[2].trim(); continue; }
    if (line.slice(0, 6) !== 'ATOM  ') continue;
    const nm = line.slice(12, 16).trim();
    const p = [parseFloat(line.slice(30, 38)), parseFloat(line.slice(38, 46)), parseFloat(line.slice(46, 54))];
    all.push({name: nm, xyz: p, res: parseInt(line.slice(22, 26), 10)});
    if (nm === 'CA') { ca.push(p); names.push(nm); }
  }
  return {meta, ca, all};
}
const centroid = P => {
  const c = [0, 0, 0];
  P.forEach(p => { c[0] += p[0]; c[1] += p[1]; c[2] += p[2]; });
  return c.map(v => v / P.length);
};
const sub = (P, c) => P.map(p => [p[0] - c[0], p[1] - c[1], p[2] - c[2]]);

/* Jacobi eigen-decomposition of a symmetric 4x4 — used for Horn's quaternion rotation. */
function jacobi4(A) {
  const a = A.map(r => r.slice()), V = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
  for (let sweep = 0; sweep < 64; sweep++) {
    let off = 0;
    for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) off += a[i][j] * a[i][j];
    if (off < 1e-22) break;
    for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) {
      if (Math.abs(a[p][q]) < 1e-24) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < 4; k++) {
        const akp = a[k][p], akq = a[k][q];
        a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
      }
      for (let k = 0; k < 4; k++) {
        const apk = a[p][k], aqk = a[q][k];
        a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
      }
      for (let k = 0; k < 4; k++) {
        const vkp = V[k][p], vkq = V[k][q];
        V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq;
      }
    }
  }
  return {vals: [a[0][0], a[1][1], a[2][2], a[3][3]], vecs: V};
}

/** Optimal PROPER rotation taking P onto Q (both already centred). Horn 1987. */
function hornRotation(P, Q) {
  let Sxx=0,Sxy=0,Sxz=0,Syx=0,Syy=0,Syz=0,Szx=0,Szy=0,Szz=0;
  for (let i = 0; i < P.length; i++) {
    const [px,py,pz] = P[i], [qx,qy,qz] = Q[i];
    Sxx+=px*qx; Sxy+=px*qy; Sxz+=px*qz;
    Syx+=py*qx; Syy+=py*qy; Syz+=py*qz;
    Szx+=pz*qx; Szy+=pz*qy; Szz+=pz*qz;
  }
  const N = [
    [Sxx+Syy+Szz, Syz-Szy,      Szx-Sxz,      Sxy-Syx],
    [Syz-Szy,     Sxx-Syy-Szz,  Sxy+Syx,      Szx+Sxz],
    [Szx-Sxz,     Sxy+Syx,     -Sxx+Syy-Szz,  Syz+Szy],
    [Sxy-Syx,     Szx+Sxz,      Syz+Szy,     -Sxx-Syy+Szz]
  ];
  const {vals, vecs} = jacobi4(N);
  let k = 0; for (let i = 1; i < 4; i++) if (vals[i] > vals[k]) k = i;
  let q = [vecs[0][k], vecs[1][k], vecs[2][k], vecs[3][k]];
  const n = Math.hypot(...q); q = q.map(v => v / n);
  const [w,x,y,z] = q;
  return [
    [w*w+x*x-y*y-z*z, 2*(x*y-w*z),     2*(x*z+w*y)],
    [2*(x*y+w*z),     w*w-x*x+y*y-z*z, 2*(y*z-w*x)],
    [2*(x*z-w*y),     2*(y*z+w*x),     w*w-x*x-y*y+z*z]
  ];
}

/** Superpose pred onto nat; returns {moved, rmsd}. */
function superpose(pred, nat) {
  const cp = centroid(pred), cn = centroid(nat);
  const P = sub(pred, cp), Q = sub(nat, cn);
  const R = hornRotation(P, Q);
  const moved = P.map(p => [
    R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2] + cn[0],
    R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2] + cn[1],
    R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2] + cn[2]
  ]);
  let s = 0;
  for (let i = 0; i < nat.length; i++) {
    s += (moved[i][0]-nat[i][0])**2 + (moved[i][1]-nat[i][1])**2 + (moved[i][2]-nat[i][2])**2;
  }
  return {moved, rmsd: Math.sqrt(s / nat.length)};
}

/* ------------------------------------------------------------------ the 3D viewer */
const VIEWERS = [];

function drawTrace(v, pts, color, radius, sphere) {
  for (let i = 0; i < pts.length - 1; i++) {
    v.addCylinder({start: {x: pts[i][0], y: pts[i][1], z: pts[i][2]},
                   end:   {x: pts[i+1][0], y: pts[i+1][1], z: pts[i+1][2]},
                   radius, color, fromCap: 2, toCap: 2});
  }
  if (sphere) pts.forEach(p => v.addSphere({center: {x: p[0], y: p[1], z: p[2]}, radius: sphere, color}));
}

/**
 * One pane: the prediction (prominent) superposed on the native (subordinate).
 * Returns the recomputed RMSD so the caller can tripwire it against the record.
 */
function mountViewer(host, tid, cfg, opts = {}) {
  const basis = opts.basis || B1;
  const natTxt = S[skey(tid, 'native', BNAT)], predTxt = S[skey(tid, cfg, basis)];
  const box = el('div', {cls: 'viewer' + (opts.small ? ' small' : ''),
                         role: 'img',
                         'aria-label': `3D structure of ${tid}, configuration ${cfg}, ${basisLabel(basis)}, predicted trace over the native`});
  host.append(box);
  if (!predTxt || !natTxt) {
    box.append(el('div', {cls: 'fallback', text: 'No structure exported for this cell — absent, not estimated.'}));
    return null;
  }
  if (!HAS3D) {
    box.append(el('div', {cls: 'fallback',
      text: '3Dmol.js (cdnjs, pinned 2.5.5) could not be loaded, so the 3D pane is unavailable. Every number on this page is computed offline and is unaffected.'}));
    return null;
  }
  const nat = parsePDB(natTxt).ca, pred = parsePDB(predTxt).ca;
  const {moved, rmsd} = superpose(pred, nat);

  const v = $3Dmol.createViewer(box, {backgroundColor: getComputedStyle(document.body).getPropertyValue('--panel').trim() || '#fff'});
  // An invisible copy of the NATIVE is what the camera frames. Every configuration of a
  // target therefore gets the identical camera, so panes are comparable by eye; framing on
  // the drawn shapes instead would rescale each pane by how far its own prediction wanders.
  const anchor = v.addModel(natTxt, 'pdb');
  v.setStyle({model: anchor}, {});
  drawTrace(v, nat, '#8a90a0', 0.10, 0.22);
  drawTrace(v, moved, '#2f5fd0', 0.26, 0.42);
  const frame = () => { v.zoomTo({model: anchor}); v.zoom(opts.small ? 0.80 : 1.15); v.render(); };
  frame();
  VIEWERS.push(v);

  const ctl = el('div', {cls: 'ctl'}, [
    el('button', {type: 'button', text: 'Reset view', title: 'Recentre and reframe on the native'})
  ]);
  ctl.firstChild.addEventListener('click', frame);
  box.append(ctl);
  return rmsd;
}

/* ------------------------------------------------------------------ chrome */
function setChrome(crumbs) {
  const c = document.getElementById('crumbs');
  c.replaceChildren();
  crumbs.forEach((b, i) => {
    if (i) c.append(el('span', {cls: 'sep', text: '/'}));
    c.append(b.href ? el('a', {href: b.href, text: b.label}) : el('span', {text: b.label}));
  });
  const bb = document.getElementById('basisBadge');
  bb.textContent = 'basis: ' + (D.meta.basis || 'unknown');
  bb.title = D.meta.basis_note || '';
  const sb = document.getElementById('statusBadge');
  if (D.meta.status && D.meta.status !== 'GENERATED') {
    sb.hidden = false;
    sb.textContent = D.meta.status + (D.meta.contains_synthetic ? ' — synthetic rows present' : '');
  } else sb.hidden = true;
  document.getElementById('footProv').textContent =
    `git ${D.meta.git_commit || '?'} · modules ${D.meta.module_hash || '?'} · target map ${D.meta.target_map_digest || '?'} · generated ${D.meta.generated || '?'}`;
}

function statCard(k, v, unit, note) {
  return el('div', {cls: 'card stat'}, [
    el('div', {cls: 'k', text: k}),
    el('div', {cls: 'v', html: v + (unit ? ` <small>${unit}</small>` : '')}),
    note ? el('div', {cls: 'n', text: note}) : null
  ]);
}

function verdictCell(v) {
  const cls = !v ? 'v-nm' : v.startsWith('BETTER') ? 'v-better' : v.startsWith('WORSE') ? 'v-worse' : 'v-nm';
  return el('td', {cls: 'verdict ' + cls, text: v || '—'});
}

/* ------------------------------------------------------------------ view: overview */
function viewOverview(main) {
  setChrome([{label: 'Overview'}]);
  const lb = D.leaderboard;
  const base = lb.find(r => r.configuration === BASELINE);
  const best = lb.filter(r => r.mean !== null)[0];

  main.append(el('h1', {text: 'Results'}));
  main.append(el('p', {cls: 'lede', html:
    `Mean full-chain C&alpha;-RMSD over ${D.meta.n_targets} cluster-disjoint dev targets. ` +
    `Every number on this page is reported on <b>two bases, always paired</b>: the ` +
    `<b>built chain</b> (${B1}) is the production result, and the <b>point cloud</b> (${B2}) ` +
    `is an explicitly labelled non-physical intermediate that reads better and is not a molecule.`}));
  main.append(el('div', {cls: 'card', style: 'margin-bottom:16px'}, [
    el('div', {cls: 'kv', html: `<b>production basis</b> ${D.meta.basis_note || ''}`}),
    el('div', {cls: 'kv', html: `<b>intermediate</b> ${D.meta.basis_secondary_note || ''}`}),
    el('div', {cls: 'kv', html: `<b>release gates</b> ${D.meta.pool_gate_rule || ''}`}),
    el('div', {cls: 'kv', html: `<b>difficulty gate</b> ${D.meta.difficulty_gate_rule || ''}`})
  ]));

  main.append(el('div', {cls: 'grid g4'}, [
    statCard('Mean RMSD — ' + pretty(best ? best.configuration : '—'),
             num(best && best.mean, 4), 'Å',
             `built chain, the production basis · ${num(best && best.mean_secondary, 4)} Å on the point cloud`),
    statCard('Baseline — ' + pretty(BASELINE), num(base && base.mean, 4), 'Å',
             'every paired comparison is against this'),
    statCard('Targets', String(D.meta.n_targets), '', 'pinned dev set, fold-clustered'),
    statCard('Configurations', String(CONFIGS.length), '', 'same targets, folds, metric, basis')
  ]));

  if (D.meta.contains_synthetic) {
    main.append(el('div', {cls: 'note', html:
      '<b>PROVISIONAL.</b> One or more configurations on this page are served by the synthetic ' +
      'pipeline-proof provider (native plus seeded displacement). Those rows are not results and ' +
      'are labelled per row. The page exists to prove the export → schema → render path end to end ' +
      'before the architecture is frozen.'}));
  }

  main.append(el('h2', {text: 'Method leaderboard'}));
  const head = ['configuration', 'n', 'mean (chain)', 'median', 'SD', 'mean (cloud)',
                'best', 'worst', `Δ vs ${BASELINE}`, 'SE', 'MDE', 'Δ/MDE', 'fold CI95',
                'W/L', 'pool gate', 'difficulty gate', 'corr', 'verdict'];
  const tb = el('tbody');
  lb.forEach(r => {
    const tr = el('tr');
    tr.append(el('td', {}, el('a', {href: `#/method/${r.configuration}`, text: pretty(r.configuration)})));
    tr.append(el('td', {cls: 'num', text: `${r.n_present}/${r.n_targets}`}));
    ['mean','median','sd','mean_secondary','best','worst'].forEach(k =>
      tr.append(el('td', {cls: 'num', text: num(r[k], 4)})));
    tr.append(el('td', {cls: 'num', text: signed(r.paired_effect)}));
    tr.append(el('td', {cls: 'num', text: num(r.se)}));
    tr.append(el('td', {cls: 'num', text: num(r.mde)}));
    tr.append(el('td', {cls: 'num', text: signed(r.effect_over_mde, 2)}));
    tr.append(el('td', {cls: 'num', text: r.ci95_fold_lo === null || r.ci95_fold_lo === undefined
      ? '—' : `[${signed(r.ci95_fold_lo, 3)}, ${signed(r.ci95_fold_hi, 3)}]`}));
    tr.append(el('td', {cls: 'num', text: r.wins === null || r.wins === undefined ? '—' : `${r.wins}/${r.losses}`}));
    tr.append(el('td', {cls: 'verdict ' + (r.pool_gate === 'PASS' ? 'v-better' : r.pool_gate === 'FAIL' ? 'v-worse' : 'v-nm'), text: r.pool_gate || '—'}));
    tr.append(el('td', {cls: 'verdict ' + (r.difficulty_gate === 'PASS' ? 'v-better' : r.difficulty_gate === 'FAIL' ? 'v-worse' : 'v-nm'), text: r.difficulty_gate || '—'}));
    tr.append(el('td', {cls: 'num', text: signed(r.corr_with_pool_best, 3)}));
    tr.append(verdictCell(r.verdict));
    tb.append(tr);
  });
  main.append(el('div', {cls: 'tablewrap'},
    el('table', {}, [el('thead', {}, el('tr', {}, head.map(h => el('th', {text: h})))), tb])));
  main.append(el('p', {cls: 'lede', style: 'margin-top:12px', html:
    'MDE = 2.8016 &times; SE, per comparison. The fold-clustered CI is shown, not the iid one: ' +
    'each of the five distogram fold models saw ~100 of the other 125 dev natives, so an iid CI ' +
    'on this instrument is anticonservative. A row whose fold CI includes zero reads NOT MEASURED.'}));

  main.append(el('h2', {text: 'Architecture'}));
  const a = D.meta.architecture || {};
  main.append(el('div', {cls: 'card kv'},
    Object.keys(a).length
      ? Object.entries(a).map(([k, v]) => el('div', {html: `<b>${k}</b> ${Array.isArray(v) ? v.join(', ') : v}`}))
      : el('div', {text: 'Architecture summary not yet supplied — the production configuration is not frozen.'})));

  main.append(el('h2', {text: 'How to check this'}));
  main.append(el('div', {cls: 'card kv'}, [
    el('div', {html: '<b>every number</b> is recomputed from the bytes of an exported PDB, not from memory'}),
    el('div', {html: '<b>every structure</b> carries PROVENANCE, MODULE_SHA, GIT_COMMIT and BASIS in its own REMARK 999 header'}),
    el('div', {html: '<b>every basis</b> is stated beside its number; the two are never mixed'}),
    el('div', {html: '<b>absent means absent</b> — a missing result is written as absent, never filled with an estimate'}),
    el('div', {html: '<b>statistics</b> come from s24/stats_lib: MDE = 2.8016 × SE per comparison, fold-clustered CI beside iid'})
  ]));

  main.append(el('h2', {text: 'Browse a target'}));
  main.append(el('p', {cls: 'lede', text:
    'Every target below opens the same-target comparison view: one target, every configuration side by side, each showing its own RMSD.'}));
  const tb2 = el('tbody');
  D.targets.forEach(t => {
    const row = byTarget.get(t.target_id) || new Map();
    const tr = el('tr');
    tr.append(el('td', {}, el('a', {href: `#/compare/${t.target_id}`, text: t.target_id})));
    tr.append(el('td', {text: t.pdb_id}));
    tr.append(el('td', {cls: 'num', text: t.length}));
    tr.append(el('td', {cls: 'num', text: t.fold}));
    tr.append(el('td', {cls: 'seq', text: t.sequence}));
    CONFIGS.forEach(c => {
      const r = row.get(c);
      tr.append(r && r.status === 'present'
        ? el('td', {cls: 'num'}, el('a', {href: `#/target/${t.target_id}/${c}`, text: num(r.rmsd, 3)}))
        : el('td', {cls: 'absent', text: 'absent'}));
    });
    tb2.append(tr);
  });
  main.append(el('div', {cls: 'tablewrap'}, el('table', {}, [
    el('thead', {}, el('tr', {}, ['target', 'pdb', 'n', 'fold', 'sequence', ...CONFIGS.map(pretty)]
      .map(h => el('th', {text: h})))), tb2])));
}

/* ------------------------------------------------------------------ view: method */
function histogram(vals, ref) {
  const W = 720, H = 150, pad = 26, hi = Math.max(8, Math.ceil(Math.max(...vals, ref || 0)));
  const nb = 24, bins = new Array(nb).fill(0);
  vals.forEach(v => bins[Math.min(nb - 1, Math.floor(v / hi * nb))]++);
  const mx = Math.max(...bins, 1);
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.setAttribute('class', 'hist');
  svg.setAttribute('role', 'img');
  svg.setAttribute('aria-label', `RMSD distribution, ${vals.length} targets, 0 to ${hi} angstrom`);
  const bw = (W - 2 * pad) / nb;
  bins.forEach((b, i) => {
    const h = (H - pad) * b / mx;
    const r = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    r.setAttribute('x', pad + i * bw + 1); r.setAttribute('y', H - pad - h);
    r.setAttribute('width', bw - 2); r.setAttribute('height', h);
    svg.append(r);
  });
  if (ref) {
    const l = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    const x = pad + (ref / hi) * (W - 2 * pad);
    l.setAttribute('x1', x); l.setAttribute('x2', x); l.setAttribute('y1', 0); l.setAttribute('y2', H - pad);
    l.setAttribute('class', 'ref'); svg.append(l);
  }
  for (let i = 0; i <= 4; i++) {
    const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    t.setAttribute('x', pad + i * (W - 2 * pad) / 4); t.setAttribute('y', H - 8);
    t.setAttribute('class', 'axis'); t.setAttribute('text-anchor', 'middle');
    t.textContent = (hi * i / 4).toFixed(1);
    svg.append(t);
  }
  return svg;
}

function viewMethod(main, cfg) {
  setChrome([{label: 'Overview', href: '#/'}, {label: pretty(cfg)}]);
  const rows = (byConfig.get(cfg) || []).slice();
  const lb = D.leaderboard.find(r => r.configuration === cfg) || {};
  const present = rows.filter(r => r.status === 'present');
  main.append(el('h1', {text: pretty(cfg)}));
  const s0 = present[0] || {};
  main.append(el('p', {cls: 'lede', text:
    `selector: ${s0.selector || '—'} · Hamiltonians: ${(s0.hamiltonians || []).join(', ') || 'none'} · ` +
    `distogram: ${s0.distogram_used === null || s0.distogram_used === undefined ? '—' : (s0.distogram_used ? 'yes' : 'no')} · ` +
    `CVaR α: ${s0.cvar_alpha === null || s0.cvar_alpha === undefined ? '—' : s0.cvar_alpha}`}));

  main.append(el('div', {cls: 'grid g2'}, [
    statCard('Mean — built chain', num(lb.mean), 'Å',
             `${lb.n_present}/${lb.n_targets} present · ${num(lb.mean_secondary)} Å point cloud (non-physical)`),
    statCard('Median', num(lb.median), 'Å', `SD ${num(lb.sd)}`),
    statCard(`Δ vs ${BASELINE}`, signed(lb.paired_effect), 'Å',
             `MDE ${num(lb.mde)} · Δ/MDE ${signed(lb.effect_over_mde, 2)}`),
    statCard('Gates', `<span style="font-size:15px">pool ${lb.pool_gate || '—'} · difficulty ${lb.difficulty_gate || '—'}</span>`,
             '', `corr(RMSD, pool_best) ${signed(lb.corr_with_pool_best, 3)} · sd ${num(lb.sd_per_target, 3)} · ${lb.verdict || ''}`),
    statCard('Verdict', `<span style="font-size:15px">${lb.verdict || '—'}</span>`, '',
             lb.ci95_fold_lo === null || lb.ci95_fold_lo === undefined ? 'no fold CI'
               : `fold CI95 [${signed(lb.ci95_fold_lo, 3)}, ${signed(lb.ci95_fold_hi, 3)}] · folds same sign ${lb.folds_same_sign}/${lb.n_folds}`)
  ]));

  if (lb.is_synthetic) main.append(el('div', {cls: 'note', html:
    '<b>SYNTHETIC.</b> This configuration is served by the pipeline-proof provider. Not a result.'}));

  main.append(el('h2', {text: 'Distribution'}));
  const base = D.leaderboard.find(r => r.configuration === BASELINE);
  main.append(el('div', {cls: 'card'}, [
    histogram(present.map(r => r.rmsd), base ? base.mean : null),
    el('div', {cls: 'legend'}, [
      el('span', {html: '<span class="sw pred"></span> targets, C&alpha;-RMSD (Å)'}),
      el('span', {html: `<span class="sw" style="background:var(--bad)"></span> ${BASELINE} mean ${num(base && base.mean, 3)} Å`})
    ])
  ]));

  main.append(el('h2', {text: 'Targets'}));
  rows.sort((a, b) => (a.rmsd === null) - (b.rmsd === null) || a.rmsd - b.rmsd);
  const tb = el('tbody');
  rows.forEach(r => {
    const b = (byTarget.get(r.target_id) || new Map()).get(BASELINE);
    const tr = el('tr');
    tr.append(el('td', {}, r.status === 'present'
      ? el('a', {href: `#/target/${r.target_id}/${cfg}`, text: r.target_id})
      : el('span', {text: r.target_id})));
    tr.append(el('td', {text: r.pdb_id}));
    tr.append(el('td', {cls: 'num', text: r.length}));
    tr.append(el('td', {cls: 'num', text: r.fold}));
    tr.append(r.status === 'present'
      ? el('td', {cls: 'num', text: num(r.rmsd, 3)})
      : el('td', {cls: 'absent', text: r.absent_reason || 'absent'}));
    tr.append(el('td', {cls: 'num', text: r.status === 'present' ? num(r.rmsd_secondary, 3) : '—'}));
    tr.append(el('td', {cls: 'num', text: b && b.status === 'present' ? num(b.rmsd, 3) : '—'}));
    tr.append(el('td', {cls: 'num', text: (r.status === 'present' && b && b.status === 'present')
      ? signed(r.rmsd - b.rmsd, 3) : '—'}));
    tr.append(el('td', {}, el('a', {href: `#/compare/${r.target_id}`, text: 'compare all'})));
    tb.append(tr);
  });
  main.append(el('div', {cls: 'tablewrap'}, el('table', {}, [
    el('thead', {}, el('tr', {}, ['target', 'pdb', 'n', 'fold', 'RMSD chain (Å)',
      'RMSD cloud (Å)', `${BASELINE} chain (Å)`, 'Δ', '']
      .map(h => el('th', {text: h})))), tb])));
}

/* ------------------------------------------------------------------ view: target */
function rmsdBar(r, recomputed, basis) {
  const shown = basis === B2 ? r.rmsd_secondary : r.rmsd;
  const other = basis === B2 ? r.rmsd : r.rmsd_secondary;
  const otherB = basis === B2 ? B1 : B2;
  const drift = (recomputed !== null && recomputed !== undefined && shown !== null)
    ? Math.abs(recomputed - shown) : 0;
  return el('div', {cls: 'rmsdbar'}, [
    el('div', {}, [
      el('div', {cls: 'lbl', text: `Cα-RMSD to native · ${basisLabel(basis)}`}),
      el('div', {cls: 'big', html: `${num(shown, 4)} <small>Å</small>`}),
      el('div', {cls: 'lbl', style: basis === B2 ? 'color:var(--bad)' : '',
        text: basis === B2
          ? `NON-PHYSICAL INTERMEDIATE — mean virtual Cα-Cα bond ${num(r.ca_bond_mean_secondary, 2)} Å against a native ${num(r.ca_bond_mean_native, 2)} Å. Not a molecule.`
          : `production result. Paired: ${num(other, 4)} Å on the ${otherB} (non-physical intermediate).`}),
      drift > 0.01 ? el('div', {cls: 'lbl', style: 'color:var(--bad)',
        text: `browser recompute disagrees by ${num(drift, 4)} Å — investigate`}) : null
    ]),
    el('div', {cls: 'meta'}, [
      ['target', `${r.target_id} · ${r.pdb_id}`],
      ['configuration', pretty(r.configuration)],
      ['basis', basis],
      ['provenance', r.provenance || '—'],
      ['length', String(r.length)],
      ['fold', String(r.fold)],
      ['selector', r.selector || '—'],
      ['Hamiltonians', (r.hamiltonians || []).join(', ') || 'none'],
      ['distogram', r.distogram_used === null || r.distogram_used === undefined ? '—' : (r.distogram_used ? 'yes' : 'no')],
      ['CVaR α', r.cvar_alpha === null || r.cvar_alpha === undefined ? '—' : String(r.cvar_alpha)],
      ['qubits', r.qubits === null || r.qubits === undefined ? '—' : String(r.qubits)],
      ['ansatz layers', r.ansatz_layers === null || r.ansatz_layers === undefined ? '—' : String(r.ansatz_layers)],
      ['candidates', r.candidate_count === null || r.candidate_count === undefined ? '—' : String(r.candidate_count)],
      ['retained', r.retained_count === null || r.retained_count === undefined ? '—' : String(r.retained_count)]
    ].map(([k, v]) => el('div', {}, [el('span', {cls: 'k', text: k}), el('span', {cls: 'v', text: v})])))
  ]);
}

function legend() {
  return el('div', {cls: 'legend'}, [
    el('span', {html: '<span class="sw pred"></span> prediction'}),
    el('span', {html: '<span class="sw nat"></span> native (subordinate)'}),
    el('span', {text: 'drag to rotate · scroll to zoom · right-drag or two-finger drag to pan'})
  ]);
}

function viewTarget(main, tid, cfg) {
  const r = (byTarget.get(tid) || new Map()).get(cfg);
  setChrome([{label: 'Overview', href: '#/'}, {label: pretty(cfg), href: `#/method/${cfg}`},
             {label: tid}]);
  if (!r) { main.append(el('div', {cls: 'note', text: 'No record for that target and configuration.'})); return; }

  let basis = B1;
  const holder = el('div');
  const paneHost = el('div');
  const toggle = el('div', {style: 'display:flex;gap:8px;margin:12px 0'});
  const draw = () => {
    holder.replaceChildren();
    paneHost.replaceChildren();
    const rec = mountViewer(paneHost, tid, cfg, {basis});
    holder.append(rmsdBar(r, rec, basis));
    [...toggle.children].forEach(b => b.setAttribute('aria-pressed', String(b.dataset.b === basis)));
  };
  [[B1, 'Built chain — production'], [B2, 'Point cloud — non-physical']].forEach(([b, lbl]) => {
    const btn = el('button', {type: 'button', text: lbl, 'aria-pressed': String(b === basis)});
    btn.dataset.b = b;
    btn.addEventListener('click', () => { basis = b; draw(); });
    toggle.append(btn);
  });
  main.append(holder, toggle, paneHost);
  draw();
  main.append(legend());

  main.append(el('h2', {text: 'Sequence'}));
  main.append(el('div', {cls: 'card seq', text: r.sequence}));

  main.append(el('h2', {text: 'This target, every configuration'}));
  const tb = el('tbody');
  CONFIGS.forEach(c => {
    const q = (byTarget.get(tid) || new Map()).get(c);
    const tr = el('tr');
    tr.append(el('td', {}, el('a', {href: `#/target/${tid}/${c}`, text: pretty(c)})));
    tr.append(q && q.status === 'present' ? el('td', {cls: 'num', text: num(q.rmsd, 4)})
                                          : el('td', {cls: 'absent', text: 'absent'}));
    tr.append(q && q.status === 'present' ? el('td', {cls: 'num', text: num(q.rmsd_secondary, 4)})
                                          : el('td', {cls: 'absent', text: '—'}));
    tb.append(tr);
  });
  main.append(el('div', {cls: 'tablewrap'}, el('table', {}, [
    el('thead', {}, el('tr', {}, [`configuration`, `${B1} (Å) — production`,
      `${B2} (Å) — non-physical`].map(h => el('th', {text: h})))), tb])));
  main.append(el('p', {style: 'margin-top:12px'},
    el('a', {href: `#/compare/${tid}`, text: 'Open the side-by-side comparison for this target →'})));

  main.append(el('h2', {text: 'Provenance'}));
  main.append(el('div', {cls: 'card kv'}, [
    [`${B1} file (production)`, r.prediction_path],
    [`${B2} file (intermediate)`, r.prediction_path_secondary],
    ['native referent', r.native_path],
    ['provenance', r.provenance], ['source', r.source], ['seed', r.seed],
    ['pool best (ORACLE)', r.pool_best === null || r.pool_best === undefined ? null : num(r.pool_best, 4) + ' Å'],
    ['mean Cα-Cα bond, chain', num(r.ca_bond_mean, 3) + ' Å'],
    ['mean Cα-Cα bond, cloud', num(r.ca_bond_mean_secondary, 3) + ' Å'],
    ['mean Cα-Cα bond, native', num(r.ca_bond_mean_native, 3) + ' Å'],
    ['git commit', r.git_commit], ['module hash', r.module_hash], ['generated', r.timestamp]
  ].map(([k, v]) => el('div', {html: `<b>${k}</b> ${v === null || v === undefined ? '—' : v}`}))));
}

/* ------------------------------------------------------------------ view: compare */
function viewCompare(main, tid) {
  const row = byTarget.get(tid) || new Map();
  const t = D.targets.find(x => x.target_id === tid) || {};
  setChrome([{label: 'Overview', href: '#/'}, {label: `${tid} — all configurations`}]);
  main.append(el('h1', {text: `${tid} · ${t.pdb_id || ''}`}));
  main.append(el('p', {cls: 'lede', html:
    `Same target, every configuration, identical camera and representation. ` +
    `length ${t.length} · fold ${t.fold}. Structures and the large number are the ` +
    `<b>built chain</b> — the production result; the small number beneath each is the ` +
    `point cloud, a non-physical intermediate. ` +
    `Every pane frames the native, so the panes are directly comparable.`}));
  main.append(el('div', {cls: 'card seq', text: t.sequence || ''}));
  main.append(legend());
  main.append(el('h2', {text: 'Configurations'}));

  const present = CONFIGS.map(c => row.get(c)).filter(r => r && r.status === 'present');
  const best = present.length ? Math.min(...present.map(r => r.rmsd)) : null;

  const grid = el('div', {cls: 'cmp'});
  main.append(grid);
  CONFIGS.forEach(c => {
    const r = row.get(c);
    const isBest = r && r.status === 'present' && best !== null && Math.abs(r.rmsd - best) < 1e-9;
    const cell = el('div', {cls: 'cell' + (isBest ? ' best' : '') + (r && r.status === 'present' ? '' : ' absent')});
    cell.append(el('h4', {}, el('a', {href: `#/target/${tid}/${c}`, text: pretty(c)})));
    cell.append(r && r.status === 'present'
      ? el('div', {cls: 'r', html: `${num(r.rmsd, 3)} <small>Å · built chain${isBest ? ' · best here' : ''}</small>`})
      : el('div', {cls: 'r', html: `<small>absent — ${(r && r.absent_reason) || 'no result'}</small>`}));
    if (r && r.status === 'present') {
      cell.append(el('div', {cls: 'lbl', style: 'font:11px var(--mono);color:var(--muted)',
        text: `${num(r.rmsd_secondary, 3)} Å point cloud (non-physical)`}));
    }
    const host = el('div');
    cell.append(host);
    grid.append(cell);
    if (r && r.status === 'present') mountViewer(host, tid, c, {small: true, basis: B1});
  });

  const nav = el('p', {style: 'margin-top:20px'});
  const idx = D.targets.findIndex(x => x.target_id === tid);
  if (idx > 0) nav.append(el('a', {href: `#/compare/${D.targets[idx - 1].target_id}`, text: '← previous target'}), el('span', {text: '   '}));
  if (idx >= 0 && idx < D.targets.length - 1) nav.append(el('a', {href: `#/compare/${D.targets[idx + 1].target_id}`, text: 'next target →'}));
  main.append(nav);
}

/* ------------------------------------------------------------------ router */
function route() {
  VIEWERS.length = 0;
  const main = document.getElementById('main');
  main.replaceChildren();
  const h = (location.hash || '#/').slice(2).split('/').filter(Boolean);
  try {
    if (h[0] === 'method' && h[1]) viewMethod(main, decodeURIComponent(h[1]));
    else if (h[0] === 'target' && h[1] && h[2]) viewTarget(main, h[1].toUpperCase(), decodeURIComponent(h[2]));
    else if (h[0] === 'compare' && h[1]) viewCompare(main, h[1].toUpperCase());
    else viewOverview(main);
  } catch (e) {
    main.append(el('div', {cls: 'note', text: 'Render error: ' + e.message}));
    throw e;
  }
  main.focus({preventScroll: true});
  window.scrollTo(0, 0);
}

if (typeof document !== 'undefined') {
  window.addEventListener('hashchange', route);
  window.addEventListener('resize', () => VIEWERS.forEach(v => { try { v.resize(); v.render(); } catch (e) {} }));
  route();
}

/* Exported under node so `test_site.js` can check the browser's superposition and RMSD
 * against the Python instrument's. The browser never uses this branch. */
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {parsePDB, superpose, hornRotation, jacobi4};
}
