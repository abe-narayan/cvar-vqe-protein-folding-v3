"""Is `core.energy`'s Legacy the real 11-term decomposition, or a proxy?

Three ways it could be a proxy and still pass a casual look:
  * the term NAMES are right but the columns are recomputed by something cheaper
  * the weights are right but applied to a different term ORDER
  * the terms are real but degenerate -- several columns carrying the same information,
    so the "11-term field" is really a one-term field wearing eleven names

So: compare the consolidated implementation column by column against the shipped
`legacy_field` + `energy_terms` on a REAL retrieval pool, check the weight vector is
applied in `TERMS` order, and measure the rank of the term matrix to show the
decomposition is not degenerate.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)


class RepShim:
    def __init__(self, PHI, PSI):
        self._phi = np.asarray(PHI, float)
        self._psi = np.asarray(PSI, float)
        self.n_states = 1
        self.n_residues = self._phi.shape[1]


def real_pool(pid="1A13", k=24):
    import numpy as np
    import peptide_db as db
    import protein_geometry as geo
    from s5.lib import B62, encode, windows_full
    from s7.amber_native import pool_for
    p = db.by_pdb(pid)
    fold = db.folds(5)[p.seq]
    _, S, PHI, PSI, _ = windows_full(pool_for(fold, p.seq), p.n)
    sim = B62[S, encode(p.seq)[None, :]].sum(1)
    idx = np.argsort(-sim, kind="stable")[:k]
    return p, PHI[idx].astype(float), PSI[idx].astype(float)


def main():
    import core
    import legacy_field as LF_ref
    import energy_terms as et
    import protein_geometry as geo

    lf = core.backend("legacy")
    out = {"legacy_backend": core.backend_name("legacy"),
           "reference_module": "legacy_field"}

    # -- term names and weights -------------------------------------------
    out["TERMS"] = list(lf.TERMS)
    out["TERMS_matches_energy_terms"] = list(lf.TERMS) == list(et.TERM_NAMES)
    out["n_terms"] = len(lf.TERMS)
    out["FITTED_WEIGHTS"] = {k: float(v) for k, v in lf.FITTED_WEIGHTS.items()}
    out["weights_match_reference"] = (
        {k: float(v) for k, v in lf.FITTED_WEIGHTS.items()}
        == {k: float(v) for k, v in LF_ref.FITTED_WEIGHTS.items()})
    out["n_distinct_weights"] = len(set(lf.FITTED_WEIGHTS.values()))
    out["weights_not_all_equal"] = len(set(lf.FITTED_WEIGHTS.values())) > 1
    out["SUBSETS_keys"] = sorted(lf.SUBSETS)

    # -- the term matrix on a real pool, consolidated vs shipped ----------
    p, PHI, PSI = real_pool()
    BB = geo.build_backbone_batch(PHI, PSI)
    T_new = np.asarray(lf.BatchLegacy(p.seq, RepShim(PHI, PSI))
                       .terms_from_coords(BB, phi=PHI, psi=PSI), float)
    T_ref = np.asarray(LF_ref.BatchLegacy(p.seq, RepShim(PHI, PSI))
                       .terms_from_coords(BB, phi=PHI, psi=PSI), float)
    out["term_matrix_shape"] = list(T_new.shape)
    out["term_matrix_max_abs_diff"] = float(np.max(np.abs(T_new - T_ref)))
    out["term_matrix_bit_identical"] = bool(np.array_equal(T_new, T_ref))

    # -- the weighted total, and that the weights go on in TERMS order ----
    wv = np.array([lf.FITTED_WEIGHTS.get(t, 0.0) for t in lf.TERMS], float)
    tot = T_new @ wv
    out["weighted_total_first5"] = [repr(float(x)) for x in tot[:5]]
    # a shuffled weight vector must give a DIFFERENT answer; if it does not, order is
    # not actually load-bearing and the "fitted weights" claim is empty
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(wv))
    while np.array_equal(wv[perm], wv):
        perm = rng.permutation(len(wv))
    out["shuffled_weights_change_total"] = bool(
        not np.allclose(T_new @ wv[perm], tot))

    # -- non-degeneracy: 11 names must carry more than 1 dimension --------
    sd = T_new.std(0)
    live = sd > 1e-12
    out["terms_with_variation"] = [t for t, l in zip(lf.TERMS, live) if l]
    out["n_terms_with_variation"] = int(live.sum())
    Z = (T_new[:, live] - T_new[:, live].mean(0)) / sd[live]
    out["term_matrix_rank"] = int(np.linalg.matrix_rank(Z, tol=1e-8))
    C = np.corrcoef(Z, rowvar=False)
    off = C[~np.eye(C.shape[0], dtype=bool)]
    out["max_abs_offdiag_correlation"] = float(np.max(np.abs(off)))
    out["decomposition_is_not_degenerate"] = bool(
        int(np.linalg.matrix_rank(Z, tol=1e-8)) >= max(2, int(live.sum()) - 1))

    print(json.dumps(out, indent=2, sort_keys=True, default=str))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "legacy_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True, default=str)
    return out

if __name__ == "__main__":
    main()
