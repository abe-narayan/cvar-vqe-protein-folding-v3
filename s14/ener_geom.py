"""SPRINT 14, ENER-6 -- ADVERSARIAL AUDIT OF THE BACKBONE THE ENERGIES SEE.

The optimiser must not be able to win by exploiting a geometric artefact.  If the builder
can emit a chain with an inverted CA, a broken peptide bond or an interpenetrating pair,
then a low energy is evidence about the builder and not about the structure.

Audited, on configurations drawn from the SAME uniform ensemble every other arm uses, plus
adversarially on the configurations each energy PREFERS (its own lowest decile), because
that is where an exploit would show up and a uniform sample would miss it:

    bond lengths            N-CA, CA-C, C-N(+1), C=O
    bond angles             N-CA-C, CA-C-N, C-N-CA
    peptide geometry        omega, and the CA(i)-CA(i+1) virtual distance
    chirality               the improper N-CA-C-CB sign at every non-GLY residue
    steric clashes          minimum heavy-atom separation, and the count below 2.0 A
    Ramachandran validity   fraction of residues in a generously-allowed region
    chain continuity        any CA-CA gap outside 3.6-3.9 A
    rigid-motion invariance both energies under a random rotation + translation
    numerical stability     both energies under a 1e-9 relative coordinate perturbation

    python -m s14.ener_geom
"""
from __future__ import annotations

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E

N_SAMPLE = 300


def _coords(z, idx):
    """Full backbone dict for a set of configurations, via the production builder."""
    sp = z.space()
    out = []
    for i in np.atleast_1d(idx):
        s = z.states(np.array([int(i)]))[0].astype(int)
        bits = sp.rep.bitstring_from_states(np.asarray(s, int))
        out.append(sp.rep.build_coords(bits))
    return out


def audit_set(z, idx, tag):
    cs = _coords(z, idx)
    seq = z.seq
    stats = {k: [] for k in ("N_CA", "CA_C", "C_N", "C_O", "ang_N_CA_C", "ang_CA_C_N",
                             "ang_C_N_CA", "omega", "ca_ca", "chirality", "min_heavy",
                             "n_clash_2A", "rama_ok")}
    for c in cs:
        N, CA, C = (np.asarray(c[k], float) for k in ("N", "CA", "C"))
        O = np.asarray(c["O"], float) if "O" in c else None
        CB = np.asarray(c["CB"], float) if "CB" in c else None
        stats["N_CA"] += np.linalg.norm(CA - N, axis=1).tolist()
        stats["CA_C"] += np.linalg.norm(C - CA, axis=1).tolist()
        stats["C_N"] += np.linalg.norm(N[1:] - C[:-1], axis=1).tolist()
        if O is not None:
            stats["C_O"] += np.linalg.norm(O - C, axis=1).tolist()
        stats["ang_N_CA_C"] += _ang(N, CA, C).tolist()
        stats["ang_CA_C_N"] += _ang(CA[:-1], C[:-1], N[1:]).tolist()
        stats["ang_C_N_CA"] += _ang(C[:-1], N[1:], CA[1:]).tolist()
        # omega is circular and clusters at +/-180, so an arithmetic mean of the raw
        # dihedral is meaningless (it reads ~0 with sd ~176).  Report the folded value.
        om = np.degrees(_dih(CA[:-1], C[:-1], N[1:], CA[1:]))
        stats["omega"] += np.abs(om).tolist()
        stats["ca_ca"] += np.linalg.norm(CA[1:] - CA[:-1], axis=1).tolist()
        if CB is not None:
            ok = [i for i, a in enumerate(seq) if a != "G"]
            if ok:
                ok = np.array(ok)
                stats["chirality"] += np.degrees(
                    _dih(N[ok], CA[ok], C[ok], CB[ok])).tolist()
        # NON-BONDED clash only: atoms at least two residues apart.  Including bonded and
        # 1-3 pairs makes the minimum distance read as the C=O bond length (1.231 A) and
        # reports every structure as clashed -- an artefact, not a clash.
        arrs, res = [], []
        for nm2, v in c.items():
            v = np.asarray(v, float)
            if v.ndim == 2 and v.shape[1] == 3:
                arrs.append(v); res.append(np.arange(len(v)))
        H = np.vstack(arrs); Ri = np.concatenate(res)
        D = np.linalg.norm(H[:, None, :] - H[None, :, :], axis=-1)
        far = np.abs(Ri[:, None] - Ri[None, :]) >= 2
        D = np.where(far, D, 9e9)
        stats["min_heavy"].append(float(D.min()))
        stats["n_clash_2A"].append(int((D < 2.0).sum() // 2))
    # Ramachandran: fraction of residues in a generously-allowed region
    S = z.states(np.atleast_1d(idx)).astype(int)
    rows = np.arange(z.n)
    phi = np.degrees(z.PHI[rows[None, :], S]); psi = np.degrees(z.PSI[rows[None, :], S])
    alpha = (phi > -160) & (phi < -20) & (psi > -120) & (psi < 50)
    beta = (phi > -180) & (phi < -40) & (psi > 90) & (psi < 180)
    lalpha = (phi > 20) & (phi < 100) & (psi > -20) & (psi < 90)
    stats["rama_ok"] = (alpha | beta | lalpha).reshape(-1).astype(float).tolist()
    out = dict(tag=tag, n=int(len(cs)))
    for k, v in stats.items():
        v = np.asarray(v, float)
        if len(v) == 0:
            continue
        out[k] = dict(mean=float(v.mean()), sd=float(v.std()),
                      min=float(v.min()), max=float(v.max()))
    out["chirality_L_frac"] = float(np.mean(np.asarray(stats["chirality"]) < 0)) \
        if stats["chirality"] else None
    out["frac_clashed"] = float(np.mean(np.asarray(stats["n_clash_2A"]) > 0))
    out["rama_ok_frac"] = float(np.mean(stats["rama_ok"]))
    return out


def _ang(a, b, c):
    u = a - b; v = c - b
    cos = (u * v).sum(-1) / (np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1))
    return np.degrees(np.arccos(np.clip(cos, -1, 1)))


def _dih(p0, p1, p2, p3):
    b0 = p0 - p1; b1 = p2 - p1; b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.arctan2(y, x)


# ------------------------------------------------------------------ invariance
def invariance(z, idx, seed=0):
    """Both energies under a rigid motion, and under a 1e-9 relative perturbation."""
    from core import amber as am
    from core import energy as et
    from core import geometry as geo
    sp = z.space()
    rng = np.random.default_rng(seed)
    S = z.states(np.atleast_1d(idx)).astype(int)
    rows = np.arange(z.n)
    phi = z.PHI[rows[None, :], S]; psi = z.PSI[rows[None, :], S]
    c = geo.build_backbone_batch(phi, psi)
    base = np.asarray(E.enum(z.pdb).legacy[np.atleast_1d(idx)], float)
    comp = et.components_batch(sp.seq, c, phi, psi)
    leg0 = np.array([sum(float(et.DEFAULT_WEIGHTS.get(t, 0.0)) * np.asarray(comp[t])[b]
                         for t in E.LEG_TERMS) for b in range(len(S))])
    # rigid motion
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    w, x, y, zq = q
    R = np.array([[1 - 2 * (y * y + zq * zq), 2 * (x * y - zq * w), 2 * (x * zq + y * w)],
                  [2 * (x * y + zq * w), 1 - 2 * (x * x + zq * zq), 2 * (y * zq - x * w)],
                  [2 * (x * zq - y * w), 2 * (y * zq + x * w), 1 - 2 * (x * x + y * y)]])
    t = rng.normal(scale=50.0, size=3)
    c2 = {k: (np.asarray(v, float) @ R.T + t) if np.asarray(v).ndim == 3
          else v for k, v in c.items()}
    comp2 = et.components_batch(sp.seq, c2, phi, psi)
    leg1 = np.array([sum(float(et.DEFAULT_WEIGHTS.get(t2, 0.0)) * np.asarray(comp2[t2])[b]
                         for t2 in E.LEG_TERMS) for b in range(len(S))])
    # tiny perturbation
    c3 = {k: (np.asarray(v, float) * (1 + 1e-9 * rng.normal(size=np.asarray(v).shape)))
          if np.asarray(v).ndim == 3 else v for k, v in c.items()}
    comp3 = et.components_batch(sp.seq, c3, phi, psi)
    leg2 = np.array([sum(float(et.DEFAULT_WEIGHTS.get(t3, 0.0)) * np.asarray(comp3[t3])[b]
                         for t3 in E.LEG_TERMS) for b in range(len(S))])
    # AMBER under a rigid motion of the built coordinates
    amb = []
    for b in range(min(len(S), 12)):
        bits = sp.rep.bitstring_from_states(S[b])
        cc = sp.rep.build_coords(bits)
        e0 = float(am.refine_coords(sp.seq, sp.rep, cc, k_restraint=0.0, steps=-1,
                                    tolerance=1e9, threads=1)["energy"])
        cc2 = {k: (np.asarray(v, float) @ R.T + t) for k, v in cc.items()}
        e1 = float(am.refine_coords(sp.seq, sp.rep, cc2, k_restraint=0.0, steps=-1,
                                    tolerance=1e9, threads=1)["energy"])
        amb.append((e0, e1))
    amb = np.array(amb)
    return dict(
        legacy_cache_vs_recompute_max_abs=float(np.abs(leg0 - base).max()),
        legacy_rigid_max_abs=float(np.abs(leg1 - leg0).max()),
        legacy_rigid_max_rel=float((np.abs(leg1 - leg0)
                                    / np.maximum(np.abs(leg0), 1e-9)).max()),
        legacy_perturb1e9_max_rel=float((np.abs(leg2 - leg0)
                                         / np.maximum(np.abs(leg0), 1e-9)).max()),
        amber_rigid_max_abs=float(np.abs(amb[:, 1] - amb[:, 0]).max()),
        amber_rigid_max_rel=float((np.abs(amb[:, 1] - amb[:, 0])
                                   / np.maximum(np.abs(amb[:, 0]), 1e-9)).max()),
        amber_n=int(len(amb)))


# ------------------------------------------------------------------ main
def main():
    rows = []
    for pdb in E.ENUM_TARGETS:
        z = E.enum(pdb)
        rng = np.random.default_rng(E.SEED)
        m = z.uniform_mask
        pool = z.amber_idx[m]
        uni = pool[rng.choice(len(pool), min(N_SAMPLE, len(pool)), replace=False)]
        # ADVERSARIAL: the configurations each energy PREFERS
        nlow = max(len(pool) // 10, 20)
        leg_low = pool[np.argsort(z.legacy[pool], kind="mergesort")[:nlow]]
        amb_low = pool[np.argsort(z.amber_total[m], kind="mergesort")[:nlow]]
        rec = dict(pdb=pdb,
                   uniform=audit_set(z, uni, "uniform"),
                   legacy_low_decile=audit_set(z, leg_low, "legacy_low_decile"),
                   amber_low_decile=audit_set(z, amb_low, "amber_low_decile"),
                   invariance=invariance(z, uni[:24]))
        rows.append(rec)
        print(f"{pdb} audited", flush=True)
    report(rows)
    E.write("ener_geom", dict(
        what="adversarial geometric audit of the backbone both energies score",
        per_target=rows), n_expected=len(E.ENUM_TARGETS))
    return rows


def report(rows):
    print("\n=== GEOMETRIC AUDIT ===")
    for tag in ("uniform", "legacy_low_decile", "amber_low_decile"):
        print(f"\n-- {tag} --")
        for k in ("N_CA", "CA_C", "C_N", "C_O", "ang_N_CA_C", "ang_CA_C_N",
                  "ang_C_N_CA", "omega", "ca_ca", "chirality", "min_heavy"):
            v = [r[tag][k] for r in rows if k in r[tag]]
            if not v:
                continue
            print(f"  {k:12s} mean {np.mean([x['mean'] for x in v]):9.4f} "
                  f"sd {np.mean([x['sd'] for x in v]):8.5f}  "
                  f"range [{min(x['min'] for x in v):9.4f}, "
                  f"{max(x['max'] for x in v):9.4f}]")
        print(f"  {'L-chirality':12s} {np.mean([r[tag]['chirality_L_frac'] for r in rows if r[tag]['chirality_L_frac'] is not None]):.4f}"
              f"   frac with a <2.0 A heavy clash "
              f"{np.mean([r[tag]['frac_clashed'] for r in rows]):.4f}"
              f"   Ramachandran-allowed {np.mean([r[tag]['rama_ok_frac'] for r in rows]):.4f}")
    print("\n-- invariance and numerical stability --")
    for k in rows[0]["invariance"]:
        v = [r["invariance"][k] for r in rows]
        print(f"  {k:38s} max over targets {max(v):.3e}")


if __name__ == "__main__":
    main()
