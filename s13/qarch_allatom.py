"""SPRINT 13 QUANTUM-ARCHITECTURE, EXPERIMENT 5 -- the SCOPE of the exact-support rule.

Section 1 established, on CA-CA distances under the ideal-geometry backbone builder, that
the support of d_ij is EXACTLY the j-i-1 residues strictly between i and j.  That rule is
now load-bearing for two other agents, so its scope has to be pinned: does it survive when
the distance is between ALL-ATOM positions, which is what AMBER's nonbonded terms actually
sum over?

MECHANISM UNDER TEST.  A residue's state is (phi_i, psi_i).  phi_i places C_i (and hence CB_i,
which is built from N_i, CA_i, C_i); psi_i places O_i and N_{i+1}.  So residue i's OWN
non-CA atoms move with residue i's torsions, while CA_i does not.  Prediction: for an
all-atom pair the support is the CLOSED interval i..j (s+1 residues), two wider than the CA
rule, and -- decisively -- s=0 (intra-residue) and s=1 (adjacent-residue) atom pairs are NOT
constant, whereas the corresponding CA-CA terms are.

If that holds, AMBER's per-term support is systematically wider and reaches down to zero
separation, which is a structural difference from Legacy's CA/CB pair terms that matters for
the Pauli weight spectrum.

Measured the same way as section 1: perturb one residue's state, take the median change of
every atom-pair distance over random base configurations, threshold at 0.10 A.

Output: `s13/results/qarch_allatom.json`.

    python -m s13.qarch_allatom
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s13 import qarch_lib as Q          # noqa: E402

ATOMS = ("N", "CA", "C", "O", "CB")
THRESH = 0.10


def support_map(pdb_id, k=8, n_base=64, seed=0, thresh=THRESH):
    sp = Q.Space(pdb_id, k)
    rng = np.random.default_rng(seed)
    n = sp.n
    base = sp.uniform(n_base, rng)
    c0, _, _ = sp.coords(base)
    # stack all atoms: index = residue * len(ATOMS) + atom
    A0 = np.stack([c0[a] for a in ATOMS], axis=2).reshape(n_base, n * len(ATOMS), 3)
    D0 = np.linalg.norm(A0[:, :, None, :] - A0[:, None, :, :], axis=-1)
    supp = np.zeros((n * len(ATOMS), n * len(ATOMS), n), bool)
    for m in range(n):
        acc = []
        for a in range(k):
            S = base.copy(); S[:, m] = a
            c, _, _ = sp.coords(S)
            A = np.stack([c[x] for x in ATOMS], axis=2).reshape(n_base, n * len(ATOMS), 3)
            D = np.linalg.norm(A[:, :, None, :] - A[:, None, :, :], axis=-1)
            acc.append(np.abs(D - D0))
        supp[:, :, m] = np.median(np.median(np.stack(acc), axis=0), axis=0) > thresh
    res = np.repeat(np.arange(n), len(ATOMS))
    atom = np.tile(np.arange(len(ATOMS)), n)
    return sp, supp, res, atom


def analyse(pdb_id, k=8, **kw):
    sp, supp, res, atom = support_map(pdb_id, k, **kw)
    n = sp.n
    P = supp.shape[0]
    iu, ju = np.triu_indices(P, 1)
    ri, rj = res[iu], res[ju]
    sep = np.abs(rj - ri)
    size = supp[iu, ju].sum(1)
    # predicted supports
    pred_ca = np.maximum(sep - 1, 0)              # the section-1 rule
    pred_all = sep + 1                            # the closed-interval rule
    ca = (atom[iu] == 1) & (atom[ju] == 1)
    noca = ~ca
    out = {"pdb": pdb_id, "n": n, "k": k, "thresh_A": THRESH, "n_atom_pairs": int(len(iu))}

    def block(mask, tag):
        if not mask.any():
            return
        out[tag] = {
            "n_pairs": int(mask.sum()),
            "mean_support": float(size[mask].mean()),
            "mean_sep": float(sep[mask].mean()),
            "frac_matching_CA_rule_s_minus_1": float((size[mask] == pred_ca[mask]).mean()),
            "frac_matching_closed_interval_s_plus_1":
                float((size[mask] == pred_all[mask]).mean()),
            "support_vs_sep": {int(s): float(size[mask & (sep == s)].mean())
                               for s in np.unique(sep[mask])},
            "frac_nonconstant_at_sep0": float((size[mask & (sep == 0)] > 0).mean())
                                        if (mask & (sep == 0)).any() else None,
            "frac_nonconstant_at_sep1": float((size[mask & (sep == 1)] > 0).mean())
                                        if (mask & (sep == 1)).any() else None,
        }
        # exact index-set check against the closed interval i..j
        exact = 0
        idx = np.flatnonzero(mask)
        for q in idx[:4000]:
            lo, hi = min(ri[q], rj[q]), max(ri[q], rj[q])
            want = np.zeros(n, bool); want[lo:hi + 1] = True
            exact += bool((supp[iu[q], ju[q]] == want).all())
        out[tag]["frac_exactly_closed_interval"] = float(exact / min(len(idx), 4000))

    block(ca, "CA_CA_pairs")
    block(noca, "pairs_involving_a_non_CA_atom")
    # per atom-type-pair detail for the shortest separations
    detail = {}
    for a in range(len(ATOMS)):
        for b in range(len(ATOMS)):
            m = (atom[iu] == a) & (atom[ju] == b) & (sep == 1)
            if m.any():
                detail[f"{ATOMS[a]}-{ATOMS[b]}"] = float(size[m].mean())
    out["mean_support_at_separation_1_by_atom_pair"] = detail
    return out


def main(pdbs=("1A13", "1CS9", "1KZ2")):
    rows = [analyse(p) for p in pdbs]
    for r in rows:
        ca = r["CA_CA_pairs"]; oa = r["pairs_involving_a_non_CA_atom"]
        print(f"  {r['pdb']}: CA-CA support {ca['mean_support']:.2f} "
              f"(s-1 rule {ca['frac_matching_CA_rule_s_minus_1']:.4f}); "
              f"non-CA support {oa['mean_support']:.2f} "
              f"(closed-interval s+1 exact {oa['frac_exactly_closed_interval']:.4f}); "
              f"sep-1 non-constant: CA {ca['frac_nonconstant_at_sep1']:.3f} "
              f"vs non-CA {oa['frac_nonconstant_at_sep1']:.3f}", flush=True)
    Q.write("qarch_allatom", {"what": "scope of the exact-support rule: CA-CA vs all-atom",
                              "rows": rows})
    return rows


if __name__ == "__main__":
    main()
