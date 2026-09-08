"""SHIFT agent -- STEP 3.  Price the chemical-shift torsion channel under the conditions
MEASURED in STEP 1, not assumed ones.

**ORACLE DIAGNOSTIC throughout.** Native torsions are corrupted / masked to synthesise a
restraint channel; nothing here is a deployable method and no number here is a headline.
What is new against `s13/tors_surface.py`:

  * the missingness mask is the **measured** TALOS-N completeness gate per target
    (s14/results/shift_avail_summary.json), not a synthetic dropout;
  * the error model is TALOS-N's **published empirical mixture** (Shen & Bax 2013, Table 1,
    34-protein validation set): 87.5 % Strong / 3.9 % Generous / 8.6 % Ambiguous, with
    "bad" rates of 2.8 % and 21 % respectively, where bad is defined by the paper as
    sqrt(dphi^2 + dpsi^2) > 60 deg, and RMSD 12.2 deg / 12.1 deg over the good predictions;
  * CONFIDENTLY WRONG torsions are priced **separately** from absent ones, which is the
    distinction a downstream VQE cares about (it can override a gap, not a confident error);
  * ambiguous (multimodal) residues are priced under four handling policies rather than
    silently collapsed to a point estimate.

    python -m s14.shift_price
"""
from __future__ import annotations

import os
import sys
import json
import time

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                      # noqa: E402
from s13 import tors_eval as EV                      # noqa: E402
from s13 import tors_surface as SF                   # noqa: E402
from s14.shift_bmrb import RESULTS                   # noqa: E402

D2R = np.pi / 180.0
SEEDS = (0, 1, 2, 3, 4)

# --- TALOS-N published empirical model (Shen & Bax 2013, J Biomol NMR 56:227, Table 1,
#     34-protein validation column).  LITERATURE-SUPPORTED, primary source read.
TN = {
    "p_strong": 0.875, "p_generous": 0.039, "p_ambiguous": 0.086,
    "bad_strong": 0.028, "bad_generous": 0.21,
    "rms_phi": 12.2, "rms_psi": 12.1,
}
TN_SIGMA = float(np.sqrt((TN["rms_phi"] ** 2 + TN["rms_psi"] ** 2) / 2.0))   # 12.15 deg

SIGMAS = (0.0, 5.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0)
COVS = (1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50)
BADS = (0.0, 0.028, 0.035, 0.05, 0.10, 0.20, 0.35)


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


# ---------------------------------------------------------------- the gross-error model
class Ramachandran:
    """Pooled empirical (phi, psi) marginal used to draw GROSS errors.  Sampling a wrong
    basin from the real Ramachandran distribution is far more faithful than a uniform draw
    on the torus, and reproduces the paper's 'bad' definition by rejection."""

    def __init__(self, nat):
        ph = np.concatenate([v[0] for v in nat.values()])
        ps = np.concatenate([v[1] for v in nat.values()])
        self.ph, self.ps = ph, ps
        self.m = len(ph)

    def draw_bad(self, rng, phi0, psi0, limit_deg=60.0):
        """Draw (phi, psi) from the marginal, conditioned on being > limit from the truth."""
        lim = limit_deg * D2R
        out_p = np.empty_like(phi0)
        out_s = np.empty_like(psi0)
        todo = np.arange(len(phi0))
        for _ in range(40):
            if len(todo) == 0:
                break
            j = rng.integers(0, self.m, len(todo))
            cp, cs = self.ph[j], self.ps[j]
            d = np.hypot(wrap(cp - phi0[todo]), wrap(cs - psi0[todo]))
            ok = d > lim
            out_p[todo[ok]] = cp[ok]
            out_s[todo[ok]] = cs[ok]
            todo = todo[~ok]
        if len(todo):                                  # fallback: rotate by 120 deg
            out_p[todo] = wrap(phi0[todo] + 120 * D2R)
            out_s[todo] = wrap(psi0[todo] + 120 * D2R)
        return out_p, out_s


# ------------------------------------------------------------------------------ masks
def measured_masks():
    with open(os.path.join(RESULTS, "shift_avail_summary.json")) as fh:
        S = json.load(fh)
    out = {}
    for r in S["rows"]:
        n = r["n"]
        m = np.array(r.get("talos_mask", [0] * n), bool) if r.get("talos_mask") \
            else np.zeros(n, bool)
        out[r["pdb"]] = (m, r["tier"])
    return out, S


def build_and_score(pdbid, phi, psi, mask):
    pf, pp, ca = EV.fill_pool(pdbid, phi, psi, mask)
    return float(I.ca_rmsd(I.build_ca(pf, pp), ca))


# ------------------------------------------------------------------------------ arms
def gauss_arm(rng, nphi, npsi, sigma):
    n = len(nphi)
    return (wrap(nphi + rng.normal(0, sigma * D2R, n)),
            wrap(npsi + rng.normal(0, sigma * D2R, n)))


def talosn_arm(rng, nphi, npsi, ram, mask, sigma=TN_SIGMA, bad_rate=None,
               amb_policy="drop", p_amb=None):
    """Return (phi, psi, keep_mask).  `mask` is the residues the completeness gate admits.

    amb_policy:
      drop        ambiguous residues carry no restraint (TALOS-N's own behaviour)
      collapse    collapse the two basins to one at random -- 50 % land in the wrong basin,
                  i.e. the naive "just take the argmax of a bimodal posterior"
      oracle2     ORACLE upper bound: the correct basin of the two is always chosen
      keep        BY CONSTRUCTION IDENTICAL to oracle2 (same code path, different seed).
                  Retained deliberately as a seed replicate, so the spread between the two
                  reported numbers is this experiment's seed-noise error bar.
    """
    n = len(nphi)
    p_amb = TN["p_ambiguous"] if p_amb is None else p_amb
    p_str = TN["p_strong"]
    p_gen = TN["p_generous"]
    tot = p_str + p_gen + p_amb
    p_str, p_gen, p_amb = p_str / tot, p_gen / tot, p_amb / tot
    u = rng.random(n)
    cls = np.where(u < p_str, 0, np.where(u < p_str + p_gen, 1, 2))   # 0 str 1 gen 2 amb

    phi = wrap(nphi + rng.normal(0, sigma * D2R, n))
    psi = wrap(npsi + rng.normal(0, sigma * D2R, n))

    bs = TN["bad_strong"] if bad_rate is None else bad_rate
    bg = TN["bad_generous"] if bad_rate is None else bad_rate
    pbad = np.where(cls == 0, bs, np.where(cls == 1, bg, bs))
    isbad = rng.random(n) < pbad
    if isbad.any():
        bp, bsx = ram.draw_bad(rng, nphi[isbad], npsi[isbad])
        phi[isbad], psi[isbad] = bp, bsx

    keep = mask.copy()
    amb = (cls == 2) & mask
    if amb.any():
        if amb_policy == "drop":
            keep = keep & ~amb
        elif amb_policy == "collapse":
            flip = rng.random(int(amb.sum())) < 0.5
            idx = np.flatnonzero(amb)[flip]
            if len(idx):
                bp, bsx = ram.draw_bad(rng, nphi[idx], npsi[idx])
                phi[idx], psi[idx] = bp, bsx
        elif amb_policy == "oracle2":
            pass                       # the near-native basin is always chosen
        elif amb_policy == "keep":
            pass
    return phi, psi, keep


# ------------------------------------------------------------------------------ driver
def main():
    t0 = time.time()
    nat = EV.native_torsions()
    ram = Ramachandran(nat)
    masks, S = measured_masks()
    tg = I.targets()
    fail18 = set(I.FAIL18)
    runnable = set(S["runnable_T3plus"])

    def agg(vals, sel=None):
        v = np.array([vals[t["pdb"]] for t in tg
                      if sel is None or t["pdb"] in sel])
        return {"mean": float(v.mean()), "median": float(np.median(v)),
                "frac_u2": float((v < 2.0).mean()), "n": int(len(v))}

    out = {"talosn_params": TN, "talosn_sigma": TN_SIGMA, "seeds": list(SEEDS)}

    # ---------------- A. sigma x coverage, MEASURED mask vs synthetic dropout
    print("A. sigma sweep on the MEASURED completeness mask (runnable T>=3, n=%d)"
          % len(runnable))
    A = {}
    for sg in SIGMAS:
        per = {}
        for t in tg:
            p = t["pdb"]
            if p not in runnable:
                continue
            np_, ps_ = nat[p]
            m = masks[p][0]
            acc = []
            for s in SEEDS:
                rng = np.random.default_rng(int(1e6 * sg + s + hash(p) % 9973))
                ph, pss = gauss_arm(rng, np_, ps_, sg)
                acc.append(build_and_score(p, ph, pss, m))
            per[p] = float(np.mean(acc))
        A["measured_sig%g" % sg] = agg(per, runnable)
        print("   sigma %5.1f  measured-mask  %.3f   <2A %.2f"
              % (sg, A["measured_sig%g" % sg]["mean"], A["measured_sig%g" % sg]["frac_u2"]))
    out["A_measured_mask"] = A

    # full coverage on the same subset, for the coverage-cost decomposition
    Afull = {}
    for sg in SIGMAS:
        per = {}
        for t in tg:
            p = t["pdb"]
            if p not in runnable:
                continue
            np_, ps_ = nat[p]
            acc = []
            for s in SEEDS:
                rng = np.random.default_rng(int(1e6 * sg + 77 + s + hash(p) % 9973))
                ph, pss = gauss_arm(rng, np_, ps_, sg)
                acc.append(float(I.ca_rmsd(I.build_ca(ph, pss), EV.pool_torsions(p)[2])))
            per[p] = float(np.mean(acc))
        Afull["full_sig%g" % sg] = agg(per, runnable)
        print("   sigma %5.1f  FULL coverage  %.3f   <2A %.2f"
              % (sg, Afull["full_sig%g" % sg]["mean"], Afull["full_sig%g" % sg]["frac_u2"]))
    out["A_full_coverage"] = Afull

    # ---------------- B. coverage sweep at sigma 12, three missingness models
    print("B. coverage sweep at sigma=%.1f, uniform / terminal / clustered (runnable set)"
          % TN_SIGMA)
    Bres = {}
    for fr in COVS:
        for mk, fn in SF.MASKS.items():
            per = {}
            for t in tg:
                p = t["pdb"]
                if p not in runnable:
                    continue
                np_, ps_ = nat[p]
                acc = []
                for s in SEEDS:
                    rng = np.random.default_rng(int(1e4 * fr * 100 + s + hash(p + mk) % 9973))
                    ph, pss = gauss_arm(rng, np_, ps_, TN_SIGMA)
                    m = fn(rng, t["n"], fr)
                    acc.append(build_and_score(p, ph, pss, m))
                per[p] = float(np.mean(acc))
            Bres["f%g_%s" % (fr, mk)] = agg(per, runnable)
        print("   cov %.2f  uniform %.3f  clustered %.3f  terminal %.3f"
              % (fr, Bres["f%g_uniform" % fr]["mean"], Bres["f%g_clustered" % fr]["mean"],
                 Bres["f%g_terminal" % fr]["mean"]))
    out["B_coverage"] = Bres

    # ---------------- C. the full TALOS-N empirical mixture, ambiguity policies
    print("C. TALOS-N empirical mixture on the measured mask")
    Cres = {}
    for pol in ("drop", "collapse", "oracle2", "keep"):
        per = {}
        for t in tg:
            p = t["pdb"]
            if p not in runnable:
                continue
            np_, ps_ = nat[p]
            m = masks[p][0]
            acc = []
            for s in SEEDS:
                rng = np.random.default_rng(int(s + hash(p + pol) % 99733))
                ph, pss, keep = talosn_arm(rng, np_, ps_, ram, m, amb_policy=pol)
                acc.append(build_and_score(p, ph, pss, keep))
            per[p] = float(np.mean(acc))
        Cres[pol] = agg(per, runnable)
        Cres[pol]["per_target"] = per
        print("   ambiguity policy %-9s  %.3f   <2A %.2f"
              % (pol, Cres[pol]["mean"], Cres[pol]["frac_u2"]))
    out["C_talosn"] = Cres

    # ---------------- D. cost of being CONFIDENTLY WRONG vs being ABSENT
    print("D. gross-error rate sweep at sigma=%.1f, full measured mask" % TN_SIGMA)
    Dres = {}
    for br in BADS:
        per = {}
        for t in tg:
            p = t["pdb"]
            if p not in runnable:
                continue
            np_, ps_ = nat[p]
            m = masks[p][0]
            acc = []
            for s in SEEDS:
                rng = np.random.default_rng(int(1e5 * br * 1000 + s + hash(p) % 9973))
                ph, pss, keep = talosn_arm(rng, np_, ps_, ram, m, bad_rate=br,
                                           amb_policy="keep", p_amb=0.0)
                acc.append(build_and_score(p, ph, pss, keep))
            per[p] = float(np.mean(acc))
        Dres["bad%g" % br] = agg(per, runnable)
        print("   bad rate %5.3f  %.3f" % (br, Dres["bad%g" % br]["mean"]))
    # matched-information control: DROP the same fraction instead of corrupting it
    for br in BADS:
        if br == 0.0:
            continue
        per = {}
        for t in tg:
            p = t["pdb"]
            if p not in runnable:
                continue
            np_, ps_ = nat[p]
            m = masks[p][0]
            acc = []
            for s in SEEDS:
                rng = np.random.default_rng(int(1e5 * br * 1000 + 31 + s + hash(p) % 9973))
                ph, pss = gauss_arm(rng, np_, ps_, TN_SIGMA)
                drop = (rng.random(t["n"]) < br) & m
                acc.append(build_and_score(p, ph, pss, m & ~drop))
            per[p] = float(np.mean(acc))
        Dres["drop%g" % br] = agg(per, runnable)
        print("   DROP rate %5.3f  %.3f  (vs corrupt %.3f)"
              % (br, Dres["drop%g" % br]["mean"], Dres["bad%g" % br]["mean"]))
    out["D_wrong_vs_absent"] = Dres

    out["what"] = ("ORACLE DIAGNOSTIC. Prices a TALOS-N-class torsion channel on the "
                   "MEASURED per-target completeness masks of s14/results/"
                   "shift_avail_summary.json, using the published TALOS-N empirical error "
                   "mixture. Native torsions are corrupted; not a method.")
    out["runnable_n"] = len(runnable)
    out["secs"] = round(time.time() - t0, 1)
    with open(os.path.join(RESULTS, "shift_price.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote s14/results/shift_price.json  [%.0fs]" % (time.time() - t0))


if __name__ == "__main__":
    main()
