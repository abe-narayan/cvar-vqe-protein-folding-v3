"""S25 / LANE Q -- Q-B. WHAT DOES THE CVaR TAIL CONTRIBUTE AT THE DEPLOYED TEMPERATURE?

PRE-REGISTERED in `s25/PREREG_Q.md` Q-B; six forks filed to and approved by the coordinator
before this ran. RULE 0 FORK LIST, each naming the alternative NOT taken:

  1 FUNCTIONAL      per-target CA-RMSD of the arm's returned structure, paired between cells.
                    NOT TAKEN: the CVaR value / free energy / state entropy as the endpoint --
                    an objective-level story is exactly what was mistaken for a contribution.
  2 BASIS           entirely INTERNAL to `s8/integrate_vqe.json`: single-candidate selection
                    (consensus medoid) out of a 128-candidate score-filtered set.
                    NOT TAKEN: comparing any number here to the 3.0483 A incumbent, which is a
                    top-75 uniform COORDINATE AVERAGE on a different operator. Never one column.
  3 READOUT         the stored `p_theta`-weighted consensus medoid, as the artefact ran it.
                    NOT TAKEN: re-running with a coordinate-average or a probability-weighted
                    CVaR readout (s23 L8 closed the latter at 4 of 4 temperatures).
  4 NORMALISATION   energies as stored: `zrank` of the shipped distogram score over the top-2^7.
                    NOT TAKEN: a raw/moment z-score (pauli-spectrum-delta-spike-artefact).
  5 NULL            the artefact's own no-circuit arms -- `argmin`, `boltz_T*`, `topfrac_*`,
                    `medoid128` -- which share the SAME candidates, energies and readout and
                    differ ONLY in the weight vector handed to `consensus_medoid`.
                    NOT TAKEN: an initialisation mean or a uniform-on-the-torus draw.
  6 THE LABEL       paired difference in mean CA-RMSD at n=126, SE, effect/MDE, iid CI beside a
                    fold-clustered CI, W/L, median beside mean, verdict from the FIXED
                    `d_harness.paired_stats` rule.
                    NOT TAKEN: labelling on marginal means (how "+0.113 A" was produced), on
                    W/L, or on a 5-cluster CI excluding zero.

DECLARED WEAKENING: I read the artefact's stored per-arm MARGINAL MEANS before filing the
forks. No paired SE, MDE, CI or verdict had been computed. Declared in the prereg and here.

FALSIFIER: if the alpha effect at T=0.3 and T=1.0 matches the T=0.1 effect in sign and size,
my prediction is wrong, the inherited "+0.113 A is the CVaR tail" attribution stands, and this
module says so.

WHY THE ENTROPY AXIS IS COMPUTABLE FOR EVERY ARM WITHOUT RE-RUNNING ANYTHING
============================================================================
`s8/integrate.py:1272-1320` hands EVERY arm the same operator: `consensus_medoid(D, o, w)`,
the member of the same 128-candidate set minimising a w-weighted mean distance. The arms
differ ONLY in w. And the energy is `E = zrank(s[o])` where `o` is the score-ORDER, so E is
the standardised ranks 1..128 -- **the same vector for every target**. Therefore:

    argmin        w = delta            ->  H = 0 bits          (exactly)
    topfrac_f     w = uniform on top k ->  H = log2(k)         (exactly)
    boltz_T       w = exp(-E/T)        ->  H computable from E (exactly, same for all targets)
    medoid128     w = uniform on 128   ->  H = 7 bits          (exactly)
    vqe_a_T       w = p_theta          ->  H stored per arm as the mean over 126 targets

So RMSD can be plotted against READOUT ENTROPY on one axis for all 18 arms, and the question
"is there one curve, or does the circuit sit off it" becomes answerable from stored data.
LIMITATION, STATED: the VQE entropies are per-arm MEANS over targets, the others are exact
constants. The curve is therefore a CELL-LEVEL relationship (n=18 cells), not a per-target one.
"""
from __future__ import annotations

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s24 import d_harness as H          # noqa: E402
from s24 import stats_lib as ST         # noqa: E402

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

V = json.load(open(os.path.join(ROOT, "s8", "integrate_vqe.json")))
PT = {k: np.asarray(v, float) for k, v in V["per_target"].items()}
FOLDS = np.asarray(V["folds"], int)
N = len(FOLDS)
ALPHAS = tuple(V["alphas"])
TEMPS = tuple(V["temps"])
DIM = int(V["dim"])

R = {"n": N, "dim": DIM, "alphas": list(ALPHAS), "temps": list(TEMPS),
     "basis": "s8 instrument: consensus-medoid SELECTION from a 128-candidate score-filtered "
              "set. NOT comparable to the 3.0483 A top-75 coordinate-average incumbent.",
     "declared_weakening": "marginal means were read before the forks were filed"}


def line(st, tag):
    return (f"    {tag:<34} {st['mean']:+.4f}  SE {st['se']:.4f}  MDE {st['mde']:.4f}  "
            f"{st['eff_over_mde']:.2f}x  {st['W']}W/{st['L']}L  med {st['median']:+.4f}  "
            f"fold {st.get('folds_same_sign', '-')}/5  {st['verdict']}")


def cmp(a_key, b_key, tag, store):
    st = H.paired_stats(PT[a_key], PT[b_key], FOLDS, seed=0, name_a=a_key, name_b=b_key)
    print(line(st, tag))
    store[tag] = st
    return st


# =============================================== 0. the inherited claim, re-verified exactly
print("\n0. THE INHERITED CLAIM, RE-VERIFIED ARITHMETICALLY")
c_a1T01 = float(PT["vqe_a1.0_T0.1"].mean())
c_a01T01 = float(PT["vqe_a0.1_T0.1"].mean())
print(f"    vqe_a1.0_T0.1  {c_a1T01:.4f}   entropy {V['entropy_bits']['vqe_a1.0_T0.1']:.4f} bits")
print(f"    vqe_a0.1_T0.1  {c_a01T01:.4f}   entropy {V['entropy_bits']['vqe_a0.1_T0.1']:.4f} bits")
print(f"    difference     {c_a1T01 - c_a01T01:+.4f}   (the inherited '+0.113 A')")
print(f"    argmin arm     {float(PT['argmin'].mean()):.4f}   "
      f"identical to vqe_a1.0_T0.1: {np.allclose(PT['argmin'], PT['vqe_a1.0_T0.1'])}")
R["inherited_claim"] = dict(
    a1_T01=c_a1T01, a01_T01=c_a01T01, difference=c_a1T01 - c_a01T01,
    entropy_a1_T01=V["entropy_bits"]["vqe_a1.0_T0.1"],
    entropy_a01_T01=V["entropy_bits"]["vqe_a0.1_T0.1"],
    a1_T01_is_argmin_per_target=bool(np.allclose(PT["argmin"], PT["vqe_a1.0_T0.1"])),
    verified=True)

# ====================================================== 1. PRIMARY: the alpha effect by T
print("\n1. PRIMARY (registered). THE ALPHA EFFECT AT EACH TEMPERATURE.")
print("   NEGATIVE means the TAIL (alpha<1) is BETTER than no tail (alpha=1).")
prim = {}
for T in TEMPS:
    print(f"  T = {T}")
    for a in (0.1, 0.25):
        cmp(f"vqe_a{a}_T{T}", f"vqe_a1.0_T{T}", f"a={a} - a=1.0  @T={T}", prim)
R["primary_alpha_effect_by_T"] = prim

print("\n   marginal means for orientation (means only, NOT the label):")
for T in TEMPS:
    row = "   ".join(f"a{a}:{float(PT[f'vqe_a{a}_T{T}'].mean()):.4f}" for a in ALPHAS)
    print(f"    T={T:<5} {row}")

# ============================================ 2. IS THE TAIL WORTH ANYTHING AS DEPLOYED?
print("\n2. THE DEPLOYED QUESTION. T = 0.3 IS THE ONLY TEMPERATURE IN VQE_LFO.")
dep = {}
cmp("vqe_a0.25_T0.3", "vqe_a1.0_T0.3", "DEPLOYED T: a=0.25 - a=1.0", dep)
cmp("vqe_a0.1_T0.3", "vqe_a1.0_T0.3", "DEPLOYED T: a=0.10 - a=1.0", dep)
R["deployed_temperature"] = dep

# the deployed table itself, and what its three alpha=1 folds cost or buy
from core.pipeline import VQE_LFO            # noqa: E402
lfo_arm = np.array([PT[f"vqe_a{VQE_LFO[int(f)][0]}_T{VQE_LFO[int(f)][1]}"][i]
                    for i, f in enumerate(FOLDS)])
print(f"\n    VQE_LFO reconstructed from core/pipeline    {lfo_arm.mean():.4f}")
print(f"    stored vqe_LFO arm                          {float(PT['vqe_LFO'].mean()):.4f}")
print(f"    reproduces per target                       "
      f"{bool(np.allclose(lfo_arm, PT['vqe_LFO']))}")
R["lfo_reconstruction"] = dict(mean=float(lfo_arm.mean()),
                               stored=float(PT["vqe_LFO"].mean()),
                               per_target_identical=bool(np.allclose(lfo_arm, PT["vqe_LFO"])))

alpha1_mask = np.isin(FOLDS, [f for f in VQE_LFO if VQE_LFO[f][0] == 1.0])
print(f"    targets on an alpha=1.0 fold (NO tail)      {int(alpha1_mask.sum())} / {N} "
      f"= {alpha1_mask.mean():.1%}")
R["share_of_targets_with_no_tail_constraint"] = float(alpha1_mask.mean())

# ORACLE counterfactual, LABELLED, not a proposal: force alpha<1 everywhere.
print("\n    ORACLE COUNTERFACTUAL -- LABELLED ORACLE, NOT A PROPOSAL, NOT A TABLE CHANGE:")
for a in (0.1, 0.25):
    st = H.paired_stats(PT[f"vqe_a{a}_T0.3"], PT["vqe_LFO"], FOLDS, seed=0)
    print(line(st, f"ORACLE force a={a},T=0.3 vs VQE_LFO"))
    R.setdefault("oracle_force_alpha", {})[f"a{a}"] = st

# ================================ 3. THE ENTROPY AXIS -- ONE CURVE, OR DOES THE CIRCUIT MOVE?
print("\n3. RMSD AGAINST REALISED READOUT ENTROPY, ALL 18 ARMS, ONE OPERATOR.")


def _z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / max(x.std(), 1e-12)


E_RANKS = _z(np.arange(1, DIM + 1, dtype=float))     # = zrank(s[o]); identical every target


def hbits(w):
    w = np.asarray(w, float)
    w = w / w.sum()
    w = w[w > 0]
    return float(-(w * np.log2(w)).sum())


ENT = {}
for a in ALPHAS:
    for T in TEMPS:
        ENT[f"vqe_a{a}_T{T}"] = float(V["entropy_bits"][f"vqe_a{a}_T{T}"])
ENT["argmin"] = 0.0
ENT[f"medoid{DIM}"] = math.log2(DIM)
for f in (0.05, 0.1, 0.25, 0.5):
    k = max(1, int(round(f * DIM)))
    ENT[f"topfrac_{f}"] = math.log2(k)
for T in TEMPS:
    ENT[f"boltz_T{T}"] = hbits(np.exp(-E_RANKS / T))

fam = {}
for k in ENT:
    fam[k] = ("VQE (circuit)" if k.startswith("vqe_") else
              "Boltzmann (exact, no circuit)" if k.startswith("boltz_") else
              "uniform top-fraction (no circuit)" if k.startswith("topfrac_") else
              "degenerate")
rows = sorted(ENT, key=lambda k: ENT[k])
print(f"    {'arm':<18} {'H_readout (bits)':>16} {'mean RMSD':>10}   family")
for k in rows:
    print(f"    {k:<18} {ENT[k]:16.4f} {float(PT[k].mean()):10.4f}   {fam[k]}")
R["entropy_curve"] = {k: dict(entropy_bits=ENT[k], mean_rmsd=float(PT[k].mean()),
                             family=fam[k]) for k in rows}

x = np.array([ENT[k] for k in rows])
y = np.array([float(PT[k].mean()) for k in rows])
# quadratic in H -- the shape the mechanism predicts (too sharp collapses, too broad dilutes)
A = np.c_[np.ones_like(x), x, x ** 2]
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
pred = A @ coef
ss = 1.0 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
print(f"\n    quadratic fit RMSD ~ H + H^2 over all {len(x)} arms:  R^2 = {ss:.4f}")
hstar = -coef[1] / (2 * coef[2]) if coef[2] != 0 else float("nan")
print(f"    turning point H* = {hstar:.3f} bits of a possible {math.log2(DIM):.0f}")
R["entropy_fit"] = dict(coef=coef.tolist(), r2=float(ss), h_star=float(hstar), n_arms=len(x))

# does the CIRCUIT sit off the curve the no-circuit arms trace?
isq = np.array([fam[k] == "VQE (circuit)" for k in rows])
coef2, *_ = np.linalg.lstsq(A[~isq], y[~isq], rcond=None)
res_q = y[isq] - A[isq] @ coef2
res_c = y[~isq] - A[~isq] @ coef2
print(f"    curve fitted on the {int((~isq).sum())} NO-CIRCUIT arms only, "
      f"then VQE arms scored against it:")
print(f"      mean residual, VQE arms        {res_q.mean():+.4f} A  "
      f"(sd {res_q.std(ddof=1):.4f}, n={int(isq.sum())})")
print(f"      residual sd of the fit itself   {res_c.std(ddof=1):.4f} A")
R["circuit_off_curve"] = dict(mean_residual_vqe=float(res_q.mean()),
                              sd_residual_vqe=float(res_q.std(ddof=1)),
                              sd_residual_nocircuit=float(res_c.std(ddof=1)),
                              n_vqe=int(isq.sum()), n_nocircuit=int((~isq).sum()))

# does ALPHA add anything once H is in the model? (the coordinator's question 2)
xa = np.array([ENT[k] for k in rows if fam[k] == "VQE (circuit)"])
ya = np.array([float(PT[k].mean()) for k in rows if fam[k] == "VQE (circuit)"])
aa = np.array([float(k.split("_T")[0][5:]) for k in rows if fam[k] == "VQE (circuit)"])
tt = np.array([float(k.split("_T")[1]) for k in rows if fam[k] == "VQE (circuit)"])
M1 = np.c_[np.ones_like(xa), xa, xa ** 2]
M2 = np.c_[M1, aa]
r1 = ya - M1 @ np.linalg.lstsq(M1, ya, rcond=None)[0]
r2 = ya - M2 @ np.linalg.lstsq(M2, ya, rcond=None)[0]
print(f"\n    over the 9 VQE cells:   RSS(H, H^2) = {(r1**2).sum():.5f}   "
      f"RSS(H, H^2, alpha) = {(r2**2).sum():.5f}")
print(f"    alpha's marginal share of the remaining variance: "
      f"{1 - (r2**2).sum() / max((r1**2).sum(), 1e-12):.3f}")
R["alpha_beyond_entropy"] = dict(rss_H=float((r1 ** 2).sum()), rss_H_alpha=float((r2 ** 2).sum()),
                                 marginal_share=float(1 - (r2 ** 2).sum() /
                                                      max((r1 ** 2).sum(), 1e-12)))
print("    correlation(H, mean RMSD) over the 9 VQE cells: "
      f"{np.corrcoef(xa, ya)[0, 1]:+.4f}")
print("    correlation(alpha, mean RMSD) over the 9 VQE cells: "
      f"{np.corrcoef(aa, ya)[0, 1]:+.4f}")
print("    correlation(T, mean RMSD)     over the 9 VQE cells: "
      f"{np.corrcoef(tt, ya)[0, 1]:+.4f}")
R["cell_correlations"] = dict(H=float(np.corrcoef(xa, ya)[0, 1]),
                              alpha=float(np.corrcoef(aa, ya)[0, 1]),
                              T=float(np.corrcoef(tt, ya)[0, 1]))

# ================================= 4. THE CIRCUIT AGAINST ITS OWN UNCONSTRAINED OPTIMUM
print("\n4. THE CIRCUIT vs ITS OWN UNCONSTRAINED OPTIMUM (the expressivity question).")
print("   At alpha = 1 the objective is  mean(E) - T H(p),  whose UNCONSTRAINED minimiser over")
print("   the simplex is EXACTLY the Boltzmann distribution exp(-E/T)/Z. So `boltz_T` is not a")
print("   loose analogy: it is the exact optimum the 21-parameter RY/CNOT state is approximating.")
exp_ = {}
for T in TEMPS:
    cmp(f"vqe_a1.0_T{T}", f"boltz_T{T}", f"circuit - exact Boltzmann @T={T}", exp_)
R["expressivity_vs_boltzmann"] = exp_

print("\n   and the whole selector against the cheapest no-circuit arms:")
noc = {}
cmp("vqe_LFO", "argmin", "VQE_LFO - argmin (shipped)", noc)
cmp("vqe_LFO", "boltz_T0.3", "VQE_LFO - Boltzmann T=0.3", noc)
cmp("vqe_LFO", "topfrac_0.5", "VQE_LFO - uniform top-64", noc)
cmp("vqe_LFO", f"medoid{DIM}", "VQE_LFO - uniform top-128", noc)
R["vs_no_circuit"] = noc

# ================================================================================ verdict
print("\n5. FALSIFIER CHECK")
e01 = prim["a=0.1 - a=1.0  @T=0.1"]["mean"]
e03 = prim["a=0.1 - a=1.0  @T=0.3"]["mean"]
e10 = prim["a=0.1 - a=1.0  @T=1.0"]["mean"]
same = (np.sign(e03) == np.sign(e01)) and (np.sign(e10) == np.sign(e01)) \
    and abs(e03) > 0.5 * abs(e01) and abs(e10) > 0.5 * abs(e01)
print(f"    alpha=0.1 vs alpha=1.0:  T=0.1 {e01:+.4f}   T=0.3 {e03:+.4f}   T=1.0 {e10:+.4f}")
print(f"    FALSIFIER FIRED AGAINST ME (effect reproduces at other T): {bool(same)}")
R["falsifier_fired_against_me"] = bool(same)

ST.save_atomic(os.path.join(OUT, "q_alpha.json"),
               dict(kind="Q-B directional re-analysis of a fixed artefact", lane="Q",
                    sprint=25, prereg="s25/PREREG_Q.md Q-B", source="s8/integrate_vqe.json",
                    results=R),
               module_file=__file__)
print("\nwrote", os.path.join(OUT, "q_alpha.json"))
