"""s26/p_ladder_test.py -- synthetic unit tests for s26/p_ladder.py.  No native RMSD is read;
every real-data check is native-free (feature equality, risk-table identity, NaN-poison).

    python s26/p_ladder_test.py
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.chdir(ROOT)

import numpy as np                                   # noqa: E402

from s26 import p_ladder as L                        # noqa: E402
from core import predict as PR                       # noqa: E402


def test_lean_trainer_matches_mlp_fit():
    """LeanTrainer on gathered rows == core.predict.MLP.fit on the materialised matrix (same seed).
    The only tolerated difference is the float64-chunked moments."""
    rng = np.random.default_rng(0)
    nres, k, dS, npairs = 300, 8, 5, 4000
    R = rng.normal(size=(nres, k)).astype(np.float32)
    S = rng.normal(size=(npairs, dS)).astype(np.float32)
    ii = rng.integers(0, nres, npairs); jj = rng.integers(0, nres, npairs)
    Y = rng.integers(0, PR.NBINS, npairs)
    X = L.LeanTrainer._rows(S, R, ii, jj, np.arange(npairs))
    d_in = X.shape[1]
    m_ref = PR.MLP(d_in, seed=3, dropout=0.0).fit(X, Y, epochs=3, batch=512)
    m_lean = L.LeanTrainer(d_in, seed=3).fit(S, R, ii, jj, Y, epochs=3, batch=512)
    assert np.allclose(m_ref.mu, m_lean.mu, atol=1e-5), np.abs(m_ref.mu - m_lean.mu).max()
    assert np.allclose(m_ref.sd, m_lean.sd, atol=1e-5), np.abs(m_ref.sd - m_lean.sd).max()
    pa, pb = m_ref.predict_proba(X), m_lean.predict_proba(X)
    err = float(np.abs(pa - pb).max())
    assert err < 1e-4, err
    print("  lean trainer == MLP.fit on synthetic data: max |dprob| %.2e" % err)


def test_features_match_shipped_construction():
    """features_for('pca32') reproduces core.predict.features(use_esm=True) column for column on
    real target sequences, from the compact tables only (never the bank)."""
    L.install_guard()
    from s12 import instrument as I
    worst = 0.0; frac = 0.0; n = 0
    for t in I.targets()[:12]:
        X, i, j = L.features_for("pca32", t["seq"], t["fold"])
        Xs, i2, j2 = PR.features(t["seq"], True)
        assert X.shape == Xs.shape, (X.shape, Xs.shape)
        assert np.array_equal(i, i2) and np.array_equal(j, j2)
        worst = max(worst, float(np.abs(X - Xs).max())); frac += float((X != Xs).mean()); n += 1
        assert L.d_in_for("pca32") == Xs.shape[1] == PR.feature_width(True)
    print("  pca32 features vs shipped construction on 12 targets: max abs %.2e, mean frac of entries differing %.4f"
          % (worst, frac / n))
    assert worst < 1e-4, worst
    Xn, _, _ = L.features_for("noesm", I.targets()[0]["seq"], 0)
    Xs0, _, _ = PR.features(I.targets()[0]["seq"], False)
    assert np.array_equal(Xn, Xs0)
    print("  noesm features == core.predict.features(use_esm=False): exact")
    for rung in ("pca128", "raw"):
        X, i, j = L.features_for(rung, I.targets()[0]["seq"], I.targets()[0]["fold"])
        assert X.shape[1] == L.d_in_for(rung), (rung, X.shape)
    print("  pca128 / raw widths: %d / %d" % (L.d_in_for("pca128"), L.d_in_for("raw")))


def test_mix_identity_and_histogram():
    from s12 import instrument as I
    t = I.targets()[0]; u = I.load_univ(t["pdb"])
    P0, i, j = L.posterior("shipped", t["pdb"], t["seq"], t["fold"])
    idx = I.pool_idx(u); Dp = I.pair_dists(np.asarray(u["W"][idx], float), i, j)
    Q, _, _ = L.posterior("mix", t["pdb"], t["seq"], t["fold"], lam=0.0, D_pool=Dp)
    r0 = L.risk_table(t["seq"], P0, i, j)["risk"]; r1 = L.risk_table(t["seq"], Q, i, j)["risk"]
    assert np.array_equal(r0, r1)
    H = L.pool_histogram(Dp)
    assert H.shape == P0.shape and np.allclose(H.sum(1), 1.0)
    dg = I.distogram(t["pdb"])
    assert np.abs(r0 - np.asarray(dg["risk"], np.float32)).max() == 0.0, "risk_table is not the shipped functional"
    print("  mix lam=0 is the incumbent bit-exactly; risk_table reproduces the cached shipped risk at 0.0")


def test_progress_definitions():
    rng = np.random.default_rng(1)
    P0 = rng.dirichlet(np.ones(17), size=20); nb = rng.integers(0, 17, 20); O = np.eye(17)[nb]
    C = np.asarray(PR.CENTRES); E0 = (P0 * C).sum(1); Dt = C[nb]
    z = L.progress(P0, P0, O, E0, E0, Dt)
    assert np.isnan(z["gam_prob"]) or abs(z["gam_prob"]) < 1e-12
    z = L.progress(O, P0, O, Dt, E0, Dt)
    assert abs(z["gam_prob"] - 1.0) < 1e-12 and abs(z["cos_prob"] - 1.0) < 1e-12
    assert abs(z["gam_loc"] - 1.0) < 1e-12 and abs(z["cos_loc"] - 1.0) < 1e-12
    Qm = 0.9 * P0 + 0.1 * O
    z = L.progress(Qm, P0, O, (Qm * C).sum(1), E0, Dt)
    assert abs(z["gam_prob"] - 0.1) < 1e-12 and abs(z["cos_prob"] - 1.0) < 1e-12
    print("  progress(): gam_eff 0 at P0, 1 at the native one-hot, 0.1 at MASS(0.1); cos 1 along the ladder")


def test_nan_poison_deployable_half():
    """score_target reads no native: poison nat_ca and rr and require bit-identical scores/cloud."""
    from s12 import instrument as I
    t = I.targets()[1]; u = I.load_univ(t["pdb"])
    P0, i, j = L.posterior("shipped", t["pdb"], t["seq"], t["fold"])
    idx = I.pool_idx(u); Wp = np.asarray(u["W"][idx], float)
    sc, o, C, _, _ = L.score_target(P0, i, j, t["seq"], Wp)
    u["nat_ca"][:] = np.nan; u["rr"][:] = np.nan
    sc2, o2, C2, _, _ = L.score_target(P0, i, j, t["seq"], np.asarray(u["W"][idx], float))
    assert np.array_equal(sc, sc2) and np.array_equal(o, o2) and np.array_equal(C, C2)
    print("  NaN-poison: the deployable half is bit-identical with the native destroyed")


def test_guard_refuses_unknown_sequence():
    L.install_guard()
    from core import data as Dd
    try:
        Dd.esm_raw("WWWWWWWWWWWWWWWWWWWWWWWWWQ")
        raise AssertionError("guard did not refuse")
    except RuntimeError:
        pass
    z = Dd.esm_raw(L.PROBE_SEQ)
    assert z[0].shape == (10, 1280) and not z[0].any()
    print("  guard: unknown sequence -> RuntimeError; width-probe sequence -> zeros")


if __name__ == "__main__":
    t0 = time.time()
    for f in (test_lean_trainer_matches_mlp_fit, test_features_match_shipped_construction,
              test_mix_identity_and_histogram, test_progress_definitions,
              test_nan_poison_deployable_half, test_guard_refuses_unknown_sequence):
        print(f.__name__); f()
    print("ALL OK (%.0fs)" % (time.time() - t0))
