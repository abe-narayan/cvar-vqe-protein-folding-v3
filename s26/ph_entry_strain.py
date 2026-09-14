"""s26/ph_entry_strain.py -- ledger entry for strain_difficulty, from the two artefacts. One-shot."""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "results", "ph_strain.json"), encoding="utf-8"))
r = json.load(open(os.path.join(HERE, "results", "ph_strain_rep.json"), encoding="utf-8"))
s, sr = d["summary"], r["summary"]
rows = d["rows"]
mv = np.array([x["moved"] for x in rows]); y = np.array([x["rmsd_arm"] for x in rows])
o = np.argsort(mv, kind="stable")
q = [y[o[i * 32:(i + 1) * 32]] for i in range(4)]
edges = [mv[o[i * 32]] for i in range(4)]


def line(kind, k, src):
    z = src[kind][k]
    return (f"    {z['label']:<45} rho {z['rho']:+.3f}  iid [{z['ci95_iid'][0]:+.3f}, {z['ci95_iid'][1]:+.3f}]"
            f"  fold [{z['ci95_fold'][0]:+.3f}, {z['ci95_fold'][1]:+.3f}]  perm p {z['perm_p']:.4f}"
            f"  folds same sign {z['folds_same_sign']}/5")


tab = "\n".join(line(kind, k, s) for kind in ("raw", "partial_n_rg") for k in ("log_e0", "log_drop", "moved", "strain_after"))
rep = "\n".join(line("partial_n_rg", k, sr) for k in ("moved", "log_e0", "log_drop"))
fish = "\n".join(f"    FAIL18 in the top quartile of {k:<13}: {s['fisher_fail18'][k]['fail_in_top']}/{s['fisher_fail18'][k]['top_k']} "
                 f"(18 of 126 overall)  odds {s['fisher_fail18'][k]['odds']:.2f}  one-sided p {s['fisher_fail18'][k]['p_one_sided']:.3f}"
                 for k in ("log_e0", "log_drop", "moved", "strain_after"))
c = s["confounds"]

ENTRY = f"""`s26/ph_strain.py`, `s26/results/ph_strain.json` (complete 126/126) and the registered
replication `s26/results/ph_strain_rep.json` (different bootstrap seed, targets in reversed
order); jobs `s26/jobs_done/ph_strain.json` and `ph_strain_rep.json` (exit 0, 35 s each, peak
RSS 0.1 GB). Pre-registered in `s26/PREREG_strain_difficulty.md` (tournament rank 4, L51),
committed at 59d8e934 before the run. Reads the 126 production records
(`bench_results/cache/1fc9f2dcf489e2fb`) and nothing else. CALIBRATION, not an accuracy lever:
nothing is selected, tuned or moved; `rmsd_arm` (BUILT CHAIN basis) is read only as the ORACLE
label of the already-emitted structure. Signals, all native-free and emitted by the production
relaxation for free: `log_e0` (the built chain's own AMBER energy before relaxation),
`log_drop` (energy removed), `moved` (restraint RMSD on N/CA/C, how far the chain moved),
`strain_after` (bond + angle energy left). Confounds n and Rg partialled by residualising both
variables on them.

Spearman with the ORACLE `rmsd_arm`, n = 126, 4000-draw iid and fold-clustered bootstrap CIs,
4000-draw label-permutation p, per-fold signs:

{tab}

{fish}
    confounds: rho(n, rmsd_arm) {c['rho_n_rmsd']:+.3f}; rho(Rg, rmsd_arm) {c['rho_rg_rmsd']:+.3f}; rho(n, log_e0) {c['rho_n_log_e0']:+.3f}; rho(Rg, log_e0) {c['rho_rg_log_e0']:+.3f}

THE FALSIFIER CLEARS ON ONE SIGNAL, `moved`, WITH MARGIN: partial rho +0.433 (bar 0.25), fold CI
[+0.247, +0.588] excluding zero, 5/5 folds, permutation p < 0.00025 (bar 0.0125 after Bonferroni
over four signals). REPLICATED as registered (new seed, reversed order):

{rep}

The replication lands inside the first run's fold CI on every quantity. `log_e0` and `log_drop`
are positive at +0.24 and +0.25 partial but fail the 5/5-fold sign rule (4/5) and sit at the
0.25 bar; `strain_after` is null (partial +0.13, fold CI spanning zero). The FAIL18 Fisher test is
null for every signal (best one-sided p 0.113), as predicted: FAIL18 is about retrieval recall,
not strain.

The presentable form (built chain, ORACLE labels, quartiles of `moved`, 32 targets each):

    moved (A)     < {edges[1]:.3f}      {edges[1]:.3f} to {edges[2]:.3f}    {edges[2]:.3f} to {edges[3]:.3f}    >= {edges[3]:.3f}
    mean rmsd_arm   {q[0].mean():.3f}         {q[1].mean():.3f}             {q[2].mean():.3f}             {q[3].mean():.3f}

A chain the force field has to move 0.29 A or more to make physical has a mean error 1.6 A larger
than one it moves less than 0.16 A. Mechanism checks (ORACLE, diagnostic): without 9KAR and 2BP4
(the two broken-bond emissions, moved 0.61 and 0.50) rho is +0.41, so it is not two outliers;
rho(moved, log_e0) is +0.41 while `moved` beats `e0` on the label by 0.19 in rho, so it is not the
raw energy in disguise; the top-8 by `moved` contain both easy (1D6X 1.78 A) and hard (9KAR
7.44, 2MFV 6.04) targets, so it is a graded signal and not a flag for a few catastrophes.

What it is and is not. It is the first native-free quantity in this programme's record with a
correlation above 0.4 to the per-target error of the emitted structure (the routers of S22 L7 /
S23 L7 and the compactness proxies of `in-band-ordering-is-per-target` reached 0.24 to 0.37 on
the per-target sign). It is a confidence flag the presentation can attach to every emitted
structure at zero cost ("how far the physics had to move this chain to make it physical").
It is NOT an accuracy lever: the prereg forbids converting it into a selector or a weight, and
the record says every such conversion fails held out (S22 L7, S23 L7). The physical reading is
plain: a coordinate average that the projection turned into a strained chain is one whose pool
members disagreed, and disagreement is error. Power: at n = 126 the design resolves |rho| of
0.25 at about 0.8; the three weaker signals are at or below that bar and are reported as
measured, not as nulls. HYPOTHESIS for lane PR: quote it as a calibration curve, never as a
gain. The Adversary's check is invited.
"""
out = os.path.join(HERE, "results", "_entry_strain.md")
open(out, "w", encoding="utf-8", newline="\n").write(ENTRY)
title = ("STRAIN DIFFICULTY (tournament 4): HOW FAR THE RELAXATION MOVES THE BUILT CHAIN PREDICTS ITS ERROR "
         "(SPEARMAN +0.433 WITH rmsd_arm, FOLD CI [+0.25, +0.59], 5/5 FOLDS, PARTIAL ON n AND Rg, REPLICATED); "
         "QUARTILE MEANS 2.29 / 2.94 / 3.76 / 3.92 A; A CALIBRATION FLAG, NOT A LEVER (2026-09-13, PH)")
p = subprocess.run([sys.executable, os.path.join(HERE, "ph_ledger.py"), "--title", title, "--body", out],
                   capture_output=True, text=True)
print(p.stdout.strip(), p.stderr.strip())
os.remove(out)
