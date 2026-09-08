"""SPRINT 13 -- ADVERSARIAL AUDIT 3: the torsion prior, the helix artefact, and a
contaminated observation pool.

C3(b) claims the 1-local torsion prior's -0.553 A advantage over random sampling is
"entirely a helix artefact", using `I.ss_of` (simplified DSSP on an ideal-geometry backbone
rebuilt from the native torsions) to split helical from non-helical targets.  Sprint 12
measured that assigner at 62% agreement with a CA-geometry assigner, so the split itself is
suspect and the correlation may be measuring the assigner.

Five attacks:

  P1  peptide_db stores phi[0] = -60 deg and psi[n-1] = -45 deg as CONSTANT PLACEHOLDERS in
      all 787 entries.  Those 1,574 fabricated angles are inside the observation pool that
      `torsion_lib2` clusters and that `coord_objval.prior_logp` counts occupancy over.
      They sit in the alpha basin.  How much of the prior's alpha mass is manufactured?
  P2  the prior is exactly 1-local, so its global minimum is the per-residue argmax.  Is
      "SA on the prior at a 5,000-evaluation budget" anything other than that deterministic
      chain?  If not, the budget-matched framing is decoration.
  P3  an INDEPENDENT helix assignment from the native CA trace alone (CA(i)-CA(i+3) and
      CA(i)-CA(i+4) distances plus the virtual CA torsion), never touching torsions or
      `core.geometry`.  Does the -0.744 correlation survive?
  P4  the confound: are helical targets simply easier for EVERY arm?  The right statistic is
      the prior's advantage over random WITHIN each stratum.
  P5  the mechanism: how close is the prior's emitted chain to a pure ideal alpha-helix?

    python -m s13.adv_prior
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
import peptide_db as pdb                   # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402
import representations as reps             # noqa: E402
from s13.coord_objval import prior_logp    # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
K = 4
OUT = {}


def _wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


def _alpha(phi, psi):
    """Alpha-basin membership: phi in [-160,-20], psi in [-120, 40] degrees."""
    p = np.degrees(_wrap(phi)); q = np.degrees(_wrap(psi))
    return (p > -160) & (p < -20) & (q > -120) & (q < 40)


# --------------------------------------------------------------------------- P1
def p1_placeholder_contamination():
    ent = pdb.load()
    rows = {"n_entries": len(ent)}
    obs, obs_clean = [], []
    for p in ent:
        n = p.n
        for i in range(n):
            obs.append((float(p.phi[i]), float(p.psi[i]), reps.residue_classes(p.seq, n)[i]))
            if i not in (0, n - 1):
                obs_clean.append((float(p.phi[i]), float(p.psi[i]), reps.residue_classes(p.seq, n)[i]))
    A = np.array([[a, b] for a, b, _ in obs]); Ac = np.array([[a, b] for a, b, _ in obs_clean])
    rows["n_obs_all"] = len(A); rows["n_obs_interior"] = len(Ac)
    rows["frac_obs_terminal"] = float(1 - len(Ac) / len(A))
    rows["alpha_frac_all"] = float(_alpha(A[:, 0], A[:, 1]).mean())
    rows["alpha_frac_interior"] = float(_alpha(Ac[:, 0], Ac[:, 1]).mean())
    ph0 = np.array([[float(p.phi[0]), float(p.psi[0])] for p in ent])
    pl = np.array([[float(p.phi[-1]), float(p.psi[-1])] for p in ent])
    rows["alpha_frac_first_residue"] = float(_alpha(ph0[:, 0], ph0[:, 1]).mean())
    rows["alpha_frac_last_residue"] = float(_alpha(pl[:, 0], pl[:, 1]).mean())
    rows["phi0_values_deg"] = sorted(set(np.round(np.degrees(ph0[:, 0]), 3).tolist()))[:4]
    rows["psi_last_values_deg"] = sorted(set(np.round(np.degrees(pl[:, 1]), 3).tolist()))[:4]

    # how much does the k=4 GENERAL table move when the placeholders are removed?
    def table(pool):
        X4 = np.column_stack([np.cos(pool[:, 0]), np.sin(pool[:, 0]),
                              np.cos(pool[:, 1]), np.sin(pool[:, 1])])
        return tl2._angles_from_centres(tl2._circ_kmeans(X4, np.ones(len(X4)), K, seed=1))
    gen = np.array([[a, b] for a, b, c in obs if c == reps.CLASS_GENERAL])
    genc = np.array([[a, b] for a, b, c in obs_clean if c == reps.CLASS_GENERAL])
    T, Tc = table(gen), table(genc)
    D = np.abs(_wrap(T[:, None, 0] - Tc[None, :, 0])) + np.abs(_wrap(T[:, None, 1] - Tc[None, :, 1]))
    rows["GENERAL_table_shift_deg_when_placeholders_dropped"] = float(np.degrees(D.min(1).max()))
    rows["GENERAL_table_all_deg"] = np.degrees(T).round(1).tolist()
    rows["GENERAL_table_interior_deg"] = np.degrees(Tc).round(1).tolist()
    # occupancy of the alpha state, with and without
    def occ(pool, tab):
        d = np.abs(_wrap(pool[:, 0][:, None] - tab[None, :, 0])) + \
            np.abs(_wrap(pool[:, 1][:, None] - tab[None, :, 1]))
        lab = d.argmin(1)
        return np.bincount(lab, minlength=len(tab)) / len(pool)
    o_all, o_int = occ(gen, T), occ(genc, Tc)
    ai = int(np.argmax(_alpha(T[:, 0], T[:, 1]) * o_all)) if _alpha(T[:, 0], T[:, 1]).any() else 0
    rows["GENERAL_occupancy_all"] = o_all.round(4).tolist()
    rows["GENERAL_occupancy_interior"] = o_int.round(4).tolist()
    rows["alpha_state_occupancy_all"] = float(o_all[_alpha(T[:, 0], T[:, 1])].sum())
    rows["alpha_state_occupancy_interior"] = float(o_int[_alpha(Tc[:, 0], Tc[:, 1])].sum())
    return rows


# --------------------------------------------------------------------------- P2/P5
def p2_prior_argmin(tg, sa_rows):
    by = {r["pdb"]: r for r in sa_rows}
    rows = []
    for t in tg:
        pdbid = t["pdb"]
        if pdbid not in by:
            continue
        seq, n = t["seq"], t["n"]
        tab = tl2.library_for(seq, K, seq)
        PHI, PSI = tab[:, :, 0], tab[:, :, 1]
        lp = prior_logp(seq, K, seq)
        u = I.load_univ(pdbid); nat = u["nat_ca"]
        idx = np.arange(n)
        s_star = lp.argmax(1)                       # the EXACT global minimiser of -sum lp
        r_star = I.ca_rmsd(I.build_ca(PHI[idx, s_star], PSI[idx, s_star]), nat)
        # P5: the pure ideal alpha helix in this library = the alpha state everywhere
        al = _alpha(PHI, PSI)
        s_al = np.where(al.any(1), al.argmax(1), lp.argmax(1))
        r_al = I.ca_rmsd(I.build_ca(PHI[idx, s_al], PSI[idx, s_al]), nat)
        # an ideal alpha helix outside the library
        r_ideal = I.ca_rmsd(I.build_ca(np.full(n, np.radians(-63.0)), np.full(n, np.radians(-42.0))), nat)
        rows.append({"pdb": pdbid, "n": n, "fold": t["fold"],
                     "prior_exact_argmin_rmsd": float(r_star),
                     "sa_prior_rmsd": float(by[pdbid]["sa_prior"]),
                     "library_alpha_chain_rmsd": float(r_al),
                     "ideal_alpha_helix_rmsd": float(r_ideal),
                     "random_legacy": float(by[pdbid]["random_legacy"]),
                     "sa_prior_equals_argmin": bool(abs(r_star - by[pdbid]["sa_prior"]) < 1e-6),
                     "alpha_states_per_residue": int(al.sum())})
    return rows


# --------------------------------------------------------------------------- P3
def ca_helicity(ca):
    """Helix fraction from the NATIVE CA TRACE ALONE.  No torsions, no core.geometry.

    Alpha-helical CA geometry: |CA_i - CA_i+3| ~ 5.0-5.6 A, |CA_i - CA_i+4| ~ 6.0-6.5 A,
    virtual CA bond angle ~ 87-93 deg and virtual CA torsion ~ +45-60 deg.  A residue is
    called H when it participates in a window satisfying the two distance windows.
    """
    ca = np.asarray(ca, float); n = len(ca)
    H = np.zeros(n, bool)
    for i in range(n - 4):
        d3a = np.linalg.norm(ca[i] - ca[i + 3]); d3b = np.linalg.norm(ca[i + 1] - ca[i + 4])
        d4 = np.linalg.norm(ca[i] - ca[i + 4])
        if 4.7 <= d3a <= 5.8 and 4.7 <= d3b <= 5.8 and 5.6 <= d4 <= 6.8:
            H[i:i + 5] = True
    return H


def ca_virtual_torsion_helicity(ca):
    """Second independent CA-only assigner: virtual torsion in [40,70] deg and angle in [80,100]."""
    ca = np.asarray(ca, float); n = len(ca)
    H = np.zeros(n, bool)
    for i in range(n - 3):
        b0 = ca[i + 1] - ca[i]; b1 = ca[i + 2] - ca[i + 1]; b2 = ca[i + 3] - ca[i + 2]
        n1 = np.cross(b0, b1); n2 = np.cross(b1, b2)
        m = np.cross(n1, b1 / np.linalg.norm(b1))
        tau = np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2)))
        ang = np.degrees(np.arccos(np.clip(np.dot(-b0, b1) / (np.linalg.norm(b0) * np.linalg.norm(b1)), -1, 1)))
        if 35.0 <= tau <= 75.0 and 78.0 <= ang <= 102.0:
            H[i:i + 4] = True
    return H


def p3_helix(tg, sa_rows):
    by = {r["pdb"]: r for r in sa_rows}
    rows = []
    for t in tg:
        pdbid = t["pdb"]
        if pdbid not in by:
            continue
        u = I.load_univ(pdbid); nat = u["nat_ca"]
        p = pdb.by_pdb(pdbid)
        ss = I.ss_of(np.asarray(p.phi, float), np.asarray(p.psi, float))   # the SPRINT's assigner
        f_ss = float(np.mean([c == "H" for c in ss]))
        f_ca = float(ca_helicity(nat).mean())
        f_vt = float(ca_virtual_torsion_helicity(nat).mean())
        rows.append({"pdb": pdbid, "fold": t["fold"], "n": t["n"],
                     "helix_ss_of": f_ss, "helix_ca_dist": f_ca, "helix_ca_torsion": f_vt,
                     "sa_prior": float(by[pdbid]["sa_prior"]),
                     "random_legacy": float(by[pdbid]["random_legacy"]),
                     "sa_legacy": float(by[pdbid]["sa_legacy"]),
                     "ORACLE_descent": float(by[pdbid]["ORACLE_descent"])})
    return rows


def _spearman(a, b):
    from scipy import stats
    r = stats.spearmanr(a, b)
    return float(r.statistic), float(r.pvalue)


def p4_strata(rows, key):
    h = np.array([r[key] for r in rows], float)
    pri = np.array([r["sa_prior"] for r in rows], float)
    ran = np.array([r["random_legacy"] for r in rows], float)
    leg = np.array([r["sa_legacy"] for r in rows], float)
    orc = np.array([r["ORACLE_descent"] for r in rows], float)
    out = {"rho_helix_vs_sa_prior": _spearman(h, pri),
           "rho_helix_vs_random": _spearman(h, ran),
           "rho_helix_vs_sa_legacy": _spearman(h, leg),
           "rho_helix_vs_ORACLE_ceiling": _spearman(h, orc),
           "rho_helix_vs_prior_advantage": _spearman(h, pri - ran)}
    for nm, sel in (("helical>0.5", h > 0.5), ("mid", (h >= 0.1) & (h <= 0.5)), ("nonhelical<0.1", h < 0.1)):
        if sel.sum() < 3:
            continue
        d = I.paired(pri[sel], ran[sel])
        out[nm] = {"n": int(sel.sum()), "sa_prior": float(pri[sel].mean()),
                   "random": float(ran[sel].mean()), "sa_legacy": float(leg[sel].mean()),
                   "ORACLE": float(orc[sel].mean()),
                   "prior_minus_random": d["mean_diff"], "ci95": d["ci95"],
                   "W": d["n_better"], "L": d["n_worse"]}
    return out


def main():
    tg = I.targets()
    sa = json.load(open(os.path.join(RESULTS, "coord_search_b5000.json")))["per_target"]
    print("P1 placeholder contamination ...")
    OUT["P1_placeholders"] = p1_placeholder_contamination()
    for k in ("frac_obs_terminal", "alpha_frac_all", "alpha_frac_interior",
              "alpha_frac_first_residue", "alpha_frac_last_residue",
              "GENERAL_table_shift_deg_when_placeholders_dropped",
              "alpha_state_occupancy_all", "alpha_state_occupancy_interior"):
        print("   ", k, OUT["P1_placeholders"][k])
    print("P2/P5 prior argmin ...")
    OUT["P2_prior_argmin"] = p2_prior_argmin(tg, sa)
    r = OUT["P2_prior_argmin"]
    print("    exact prior argmin mean RMSD  ", np.mean([x["prior_exact_argmin_rmsd"] for x in r]).round(4))
    print("    SA-on-prior mean RMSD         ", np.mean([x["sa_prior_rmsd"] for x in r]).round(4))
    print("    identical on                  ", sum(x["sa_prior_equals_argmin"] for x in r), "/", len(r))
    print("    library alpha-chain mean RMSD ", np.mean([x["library_alpha_chain_rmsd"] for x in r]).round(4))
    print("    ideal alpha-helix mean RMSD   ", np.mean([x["ideal_alpha_helix_rmsd"] for x in r]).round(4))
    print("P3 independent helix assignment ...")
    OUT["P3_helix_rows"] = p3_helix(tg, sa)
    hr = OUT["P3_helix_rows"]
    a = np.array([x["helix_ss_of"] for x in hr]); b = np.array([x["helix_ca_dist"] for x in hr])
    c = np.array([x["helix_ca_torsion"] for x in hr])
    OUT["P3_assigner_agreement"] = {
        "mean_helix_ss_of": float(a.mean()), "mean_helix_ca_dist": float(b.mean()),
        "mean_helix_ca_torsion": float(c.mean()),
        "pearson_ss_vs_cadist": float(np.corrcoef(a, b)[0, 1]),
        "pearson_ss_vs_catorsion": float(np.corrcoef(a, c)[0, 1]),
        "binary_agreement_at_0.5": float(((a > 0.5) == (b > 0.5)).mean())}
    print("   ", OUT["P3_assigner_agreement"])
    OUT["P4_strata"] = {key: p4_strata(hr, key) for key in
                        ("helix_ss_of", "helix_ca_dist", "helix_ca_torsion")}
    for key, v in OUT["P4_strata"].items():
        print(f"\n  --- stratified by {key}")
        for kk, vv in v.items():
            print("   ", kk, vv)
    json.dump(OUT, open(os.path.join(RESULTS, "adv_prior.json"), "w"), indent=1)
    print("\nwrote adv_prior.json")


if __name__ == "__main__":
    main()
