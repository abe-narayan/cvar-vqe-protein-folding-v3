"""s26/ph_entry_validity.py -- ledger entry for the heavy-atom validity axis. One-shot."""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
log = open(os.path.join(HERE, "logs", "ph_validity_report.log"), encoding="utf-8").read()


def block(label):
    i = log.find("  " + label)
    assert i >= 0, label
    j = log.find("VERDICT", i)
    j = log.find("\n", j)
    return log[i:j]


WANT = ["n_clash_2A: relaxed minus built (negative = lower after)",
        "n_clash_2p6A: relaxed minus built (negative = lower after)",
        "min_heavy: relaxed minus built (negative = lower after)",
        "bond_strain: relaxed minus built (negative = lower after)",
        "angle_strain: relaxed minus built (negative = lower after)",
        "omega_dev: relaxed minus built (negative = lower after)",
        "rama_favoured: relaxed minus built (negative = lower after)"]
blocks = "\n".join("    " + block(w).replace("\n", "\n    ").rstrip() for w in WANT)
s = json.load(open(os.path.join(HERE, "results", "ph_validity.json"), encoding="utf-8"))["summary"]
AX = ("n_clash_2A", "n_clash_2p6A", "min_heavy", "bond_strain", "angle_strain", "rama_favoured", "rama_outlier",
      "cis_frac", "chirality_L_frac", "omega_dev")
tab = "\n".join(f"    {ax:<17} built {s[ax]['before']['mean']:8.4f}   relaxed {s[ax]['after']['mean']:8.4f}   helix {s[ax]['helix']['mean']:8.4f}"
                for ax in AX)

ENTRY = f"""`s26/ph_validity.py run | report`, `s26/results/ph_validity.json` (complete 126/126), job
`s26/jobs_done/ph_validity.json` (exit 0, 2568 s, peak RSS 0.274 GB, AMBER; held 620 s at the
cap). Pre-registered in `s26/PREREG_validity_axis.md` before the run. Native-free: the panel
(`s16.energy_lib.panel`) reads no native; folds enter only the CI.

GATE, per target: the production relaxation (`refine_coords(k=10, steps=0)` on the built chain
from the cached `phi`, `psi`, exactly `core.pipeline._relax_inner`) re-run here reproduces the
cached `amber_ca` and `amber_e1` on 126 of 126 targets to max |dCA| = 0.0 A and max |dE| = 0.0
kcal/mol. The panel below therefore describes the DEPLOYED emission before and after its own
relaxation, not a re-implementation. Zero-information reference: the constant alpha-helix
(phi -63, psi -42; rama 1.000 and zero clashes by construction, S16 L27), printed beside both.

    axis              built chain   relaxed chain   constant helix
{tab}
    targets with any heavy-atom pair below 2.0 A: built 34, relaxed 1;  relaxed bond_strain above 0.05 on 2 targets (2BP4, 9KAR)

`ST.fmt` verbatim, relaxed minus built (negative = lower after), the seven axes that move:

{blocks}

READING, the numbers "validity step" rests on. (1) The relaxation removes the builder's
clashes: heavy-atom pairs below 2.0 A fall from 0.444 per target (34 targets affected) to
0.008 (1 target), fold CI [-0.555, -0.317], 34W/0L/92T; contacts below 2.6 A fall from 3.45
to 0.12 per target, 81W/0L/45T; the closest heavy-atom pair moves from 2.36 to 2.78 A, 115 of
126 targets. The registered falsifier ("clashes not halved") does not fire. (2) It pays in
covalent geometry, which the built chain had ideal by construction: bond strain 0.0 to 1.3%
relative (0/126 unchanged; 12% on 2BP4), angle strain 0.0 to 2.5%, omega non-planarity 0.0 to
6.6 degrees mean (48 degrees on 1D6X; the cis fraction rises to 0.3%, two targets). (3) It buys
no Ramachandran: favoured 0.930 to 0.911, outliers 0.024 to 0.033, both NOT MEASURED. This is
the S16 L27 finding on the production input: the restraint k = 10 on N/CA/C holds the backbone
torsions where they were and lets the force field fix the packing by bending bonds and
angles. (4) The constant helix beats the relaxed chain on every axis (0 clashes, 3.08 A minimum
separation, ideal covalent geometry, rama 1.000): a validity statistic a zero-information
reference maximises is not evidence about the force field, so the defensible statement is
the conjunction, as S16 wrote it: "removes the builder's clashes (34 targets to 1) while moving
the CA trace 0.220 A and holding the torsions, at a covalent price of 1.3% bond and 2.5% angle
strain and 6.6 degrees of omega".

THE SENTENCE, for `s26/C3_RESULT.md` and the presentation: the relaxation is a validity step
that turns 34 emissions with a sub-2 A heavy-atom overlap into 1, and 125 of 126 energies above
the 1000 kcal/mol gate into converged ones, at the price of 0.021 A of accuracy, 1.3% bond and
2.5% angle strain, 6.6 degrees of peptide-bond non-planarity, and a broken virtual bond on 2BP4
and 9KAR. POWER: n = 126, SEs 0.0005 to 0.36 on the axes; every claimed change is above 2x its
MDE with 5/5 folds; the three nulls (rama favoured 0.84x, outliers 0.55x, cis 0.46x) are
underpowered for effects below their MDEs (0.023, 0.015, 0.007) and measured against anything
larger. Descriptive; no replication is due.
"""
out = os.path.join(HERE, "results", "_entry_validity.md")
open(out, "w", encoding="utf-8", newline="\n").write(ENTRY)
title = ("THE VALIDITY AXIS OF THE PRODUCTION RELAXATION, MEASURED: HEAVY-ATOM CLASHES 0.44 PER TARGET (34 TARGETS) "
         "TO 0.008 (1), CLOSEST PAIR 2.36 TO 2.78 A, AT 1.3% BOND / 2.5% ANGLE STRAIN AND 6.6 DEG OF OMEGA, NO "
         "RAMACHANDRAN GAIN; RE-RUN BIT-IDENTICAL TO THE CACHE ON 126/126 (2026-09-14, PH)")
p = subprocess.run([sys.executable, os.path.join(HERE, "ph_ledger.py"), "--title", title, "--body", out],
                   capture_output=True, text=True)
print(p.stdout.strip(), p.stderr.strip())
os.remove(out)
