"""tests/test_s29_D.py -- lane D (the Adversary): the cost-RMSD meter `s29/s29_D_cost_audit.py`.

1. The S28 values for the SHIPPED cost reproduce to the third decimal THROUGH THE METER'S OWN
   CODE PATH: ladder rho -0.182 (CA) / -0.402 (built chain, from C2's stored rows), gradient
   cosine -0.034, the native's pool percentile 0.369 (surrogate S~; 0.368 under the shipped
   lookup), pref(ORACLE circ_best vs PROD) 0.206 (CA) / 0.071 (chain) with the pool-member
   control's pct(PROD) 0.126 / 0.020 and contrasts +0.080 / +0.051 (S28-L23b, S28-L30,
   S28-L35, S28-L36, S28-L48; artefacts under `s27/results/`).  Skipped when the S28 artefacts
   or the ladder cache are absent.
2. The ladder cache reproduces the S28 structures: PROD is the deployed average, circ_best and
   NATIVE are lane A's, and every rung's ORACLE RMSD equals C2's stored `oracle_rmsd` (the
   regenerated sub0 / circ_s0 / controls are the same structures C2 scored).
3. The callable path: a native-free callable equal to S~ meters identically to the named
   DIS_SURR (values, percentile, preferences) and its finite-difference cosine is within 1e-3
   of the analytic one.
4. NaN-poison: a callable that reads ctx.nat_ca gets NaN and the meter REFUSES it.
5. Sign conventions on a synthetic target (monkeypatched loaders): a cost equal to the ORACLE
   RMSD meters at rho +1 on every ladder, cosine +1, pref 1.0 on every better rung, native
   percentile 0; the negated cost meters at rho -1, cosine -1, pref 0.  A meter whose sign
   were wrong would pass a bad objective as good; this pins the conventions.

Synthetic data for 5; the frozen S28 artefacts (read only) for 1 to 4.
"""
from __future__ import annotations

import json
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I            # noqa: E402
from s27 import s28_A_amp as A             # noqa: E402
from s27 import s28_A2_local as A2         # noqa: E402
from s29 import s29_D_cost_audit as M      # noqa: E402

HAVE_S28 = (os.path.exists(os.path.join(A.STRUCTS, "1A13.npz")) and os.path.exists(M.C2.CA_ROWS)
            and os.path.exists(M.C2.CHAIN_ROWS) and os.path.exists(os.path.join(M.S27_CACHE, "1A13.npz")))
HAVE_CACHE = HAVE_S28 and all(os.path.exists(M.cache_path(p)) for p in ("1A13", "2L7T", "9WXK"))
needs_s28 = pytest.mark.skipif(not HAVE_S28, reason="S28 artefacts absent")
needs_cache = pytest.mark.skipif(not HAVE_CACHE, reason="ladder cache absent: run `python s29/s29_D_cost_audit.py build-cache`")
PROBE = ("1A13", "2L7T", "9WXK")


def _full_cache():
    from s25 import phys_lib as P
    return all(os.path.exists(M.cache_path(p)) for p in P.targets())


# ------------------------------------------------------------------ 1. the S28 anchors
@needs_cache
def test_shipped_cost_reproduces_S28_values_at_CA_level():
    if not _full_cache():
        pytest.skip("ladder cache incomplete")
    s = M.run_meter("DIS", basis="ca", quiet=True, save=False)
    assert s["n"] == 126
    assert abs(s["ladder_rho"]["S28"]["mean"] - (-0.182)) < 5e-4, s["ladder_rho"]["S28"]["mean"]        # S28-L35 (C2 table, DIS ladder rho)
    assert abs(s["cosine"]["mean"] - (-0.034)) < 5e-4, s["cosine"]["mean"]                              # S28-L23b
    assert abs(s["cosine"]["random_ref_mean_abs"] - 0.140) < 5e-4
    assert abs(s["cosine"]["spearman_cos_vs_prod_rmsd"] - (-0.372)) < 5e-4
    assert abs(s["pref"]["circ_best"]["mean"] - 0.206) < 5e-4                                           # S28-L35
    assert abs(s["pref"]["NATIVE"]["mean"] - 0.222) < 5e-4
    m = s["pool_member_control"]
    assert abs(m["pref_pool_member_vs_prod"] - 0.126) < 5e-4                                            # S28-L36 pct(PROD)
    assert abs(m["contrast"]["effect"] - 0.080) < 5e-4                                                  # S28-L36 / L37: +0.0799, 0.83x
    assert abs(m["contrast"]["effect_over_mde"] - 0.83) < 5e-3
    assert abs(m["h2h_beats_pool_member"]["circ_best"] - 0.626) < 5e-4
    assert abs(m["h2h_beats_pool_member"]["NATIVE"] - 0.632) < 5e-4
    assert abs(s["native_pctile"]["mean"] - 0.368) < 5e-4                                                # shipped lookup (this meter)
    # the surrogate S~ gives S28-L30's 36.9th percentile, 99/126 and 6/126
    s2 = M.run_meter("DIS_SURR", basis="ca", quiet=True, save=False)
    assert abs(s2["native_pctile"]["mean"] - 0.369) < 5e-4, s2["native_pctile"]["mean"]
    assert s2["native_pctile"]["n_prod_below_native"] == 99 and s2["native_pctile"]["n_prod_below_pool_min"] == 6
    assert abs(s2["native_pctile"]["mean_f_native_posed"] - 2.089) < 5e-4
    assert abs(s2["cosine"]["mean"] - (-0.034)) < 5e-4


@needs_s28
def test_shipped_cost_reproduces_S28_values_on_the_built_chain_from_C2_rows():
    s = M.run_meter("DIS", basis="chain-s28rows", quiet=True, save=False)
    assert s["n"] == 126
    assert abs(s["ladder_rho"]["S28"]["mean"] - (-0.402)) < 5e-4, s["ladder_rho"]["S28"]["mean"]        # S28-L48, REPORT_S28 s.12
    assert abs(s["pref"]["circ_best"]["mean"] - 0.071) < 5e-4                                           # S28-L48 DIS@chain
    m = s["pool_member_control"]
    assert abs(m["pref_pool_member_vs_prod"] - 0.020) < 5e-4
    assert abs(m["contrast"]["effect"] - 0.051) < 5e-4                                                  # S28-L48 / L49: +0.0513, 0.80x
    assert abs(m["contrast"]["effect_over_mde"] - 0.80) < 5e-3
    assert abs(m["h2h_beats_pool_member"]["circ_best"] - 0.632) < 5e-4
    assert s["cosine"] is None and s["native_pctile"] is None                                            # not defined on this basis


# ------------------------------------------------------------------ 2. the cache is the S28 ladder
@needs_cache
def test_ladder_cache_reproduces_the_S28_structures():
    from s24 import d_harness as H
    for pdb in PROBE:
        cand, dis, top, dg, frame, sur = M.load_target(pdb)
        r = M.load_rungs(pdb)
        z = np.load(os.path.join(A.STRUCTS, f"{pdb}.npz"))
        Cd, _ = H.readout_uniform(cand, top)
        assert np.abs(r["S"]["PROD"] - Cd).max() < 1e-10
        assert np.abs(r["S"]["circ_best"] - z["oracle_circ"]).max() == 0.0
        assert np.abs(r["S"]["NATIVE"] - z["oracle_aff500"]).max() == 0.0
        crow = M._c2_row(pdb, "ca")
        for k in r["oracle_rmsd"]:
            if k in crow["oracle_rmsd"]:
                assert abs(r["oracle_rmsd"][k] - crow["oracle_rmsd"][k]) < 1e-6, (pdb, k)
            assert abs(r["oracle_rmsd"][k] - float(I.ca_rmsd(r["S"][k], cand.nat_ca))) < 1e-9
        assert r["meta"]["dev_circ_s0"] < 1e-6 and r["meta"]["dev_sub0"] < 1e-6
        assert max(r["meta"]["dev_vs_c2_dis"].values()) < 1e-5


# ------------------------------------------------------------------ 3. the callable path
def _surr_cost(W, ctx):
    sur = A.Surrogate(ctx.dg, ctx.n)
    return np.array([sur.value_grad(w)[0] for w in np.asarray(W, float)])


def _native_reading_cost(W, ctx):
    return np.array([float(I.ca_rmsd(w, ctx.nat_ca)) for w in np.asarray(W, float)])


@needs_cache
def test_callable_path_matches_the_named_surrogate_and_fd_gradient_is_accurate():
    sys.modules[__name__] = sys.modules[__name__]
    c = M.Cost(__name__ + ":_surr_cost")
    a = M.Cost("DIS_SURR")
    assert c.kind == "callable" and not c.analytic and a.analytic
    for pdb in PROBE:
        r1 = M.meter_target(pdb, c); r2 = M.meter_target(pdb, a)
        for k in M.RUNGS:
            assert abs(r1["f"][k] - r2["f"][k]) < 1e-12
        assert abs(r1["native_pctile"] - r2["native_pctile"]) < 1e-12
        assert r1["pref"] == r2["pref"]
        assert abs(r1["cos"] - r2["cos"]) < 1e-3, (pdb, r1["cos"], r2["cos"])
        assert r1["grad_method"].startswith("central finite differences") and r2["grad_method"].startswith("analytic")
        assert not np.isfinite(r1["rho"]["S28"]) or abs(r1["rho"]["S28"] - r2["rho"]["S28"]) < 1e-12


# ------------------------------------------------------------------ 4. NaN poison
@needs_cache
def test_meter_refuses_a_cost_that_reads_the_native():
    c = M.Cost(__name__ + ":_native_reading_cost")
    r = M.meter_target("1A13", c, need_grad=False)
    assert r["nan_on_rungs"] == list(M.RUNGS)            # every rung NaN: ctx.nat_ca is poisoned
    with pytest.raises(ValueError, match="POISONED"):
        M.run_meter(__name__ + ":_native_reading_cost", pdbs=["1A13"], quiet=True, save=False)


# ------------------------------------------------------------------ 5. sign conventions, synthetic
def _synthetic(monkeypatch, sign=+1.0):
    rng = np.random.default_rng(3)
    n, k = 9, 60
    nat = np.cumsum(rng.normal(size=(n, 3)) * 3.0, 0)
    prod = nat + rng.normal(scale=1.2, size=(n, 3))
    S = {}
    scales = {"circ_opt": 1.5, "RAND_SIGNED[0]": 2.5, "PROD": None, "GAUSS_MATCHED[0]": 2.2, "GAUSS_0.3[0]": 1.0,
              "sub0": 0.6, "circ_s0": 0.4, "circ_best": 0.25, "NATIVE": 0.0}
    for kk, sc in scales.items():
        S[kk] = prod if sc is None else nat + rng.normal(scale=sc, size=(n, 3))
    S["NATIVE"] = nat.copy()
    W = nat[None] + rng.normal(scale=2.0, size=(k, n, 3))
    orr = {kk: float(I.ca_rmsd(S[kk], nat)) for kk in S}
    cand = SimpleNamespace(pdb="SYN", seq="A" * n, n=n, fold=0, k=k, W=W, PHI=np.zeros((k, n)), PSI=np.zeros((k, n)),
                           nat_ca=nat, oracle_rr=np.linalg.norm(W - nat[None], axis=(1, 2)))
    frame = SimpleNamespace(Wp=W, ref=W[0], k=k, n=n, top=np.arange(10))
    i, j = I.pair_index(n)
    grid = np.arange(2.0, 22.0, 0.05)
    dg = {"risk": np.zeros((len(i), len(grid))), "grid": grid, "i": i, "j": j}
    monkeypatch.setattr(M, "load_target", lambda pdb: (cand, np.zeros(k), np.arange(10), dg, frame, None))
    monkeypatch.setattr(M, "load_rungs", lambda pdb, need_chain=False: {"S": S, "oracle_rmsd": orr, "meta": {}})
    monkeypatch.setattr(M.C2, "RAMA", np.zeros((1, 20, 36, 36)))

    def oracle_cost(Wc, ctx):                     # reads the TEST's native, not ctx (which is poisoned)
        return sign * np.array([float(I.ca_rmsd(w, nat)) for w in np.asarray(Wc, float)])
    return oracle_cost, orr


def test_sign_conventions_on_a_synthetic_target(monkeypatch):
    f, orr = _synthetic(monkeypatch, +1.0)
    monkeypatch.setattr(sys.modules[__name__], "_syn_cost", f, raising=False)
    c = M.Cost(__name__ + ":_syn_cost")
    r = M.meter_target("SYN", c)
    for nm in M.LADDERS:
        assert abs(r["rho"][nm] - 1.0) < 1e-12, (nm, r["rho"][nm])
    assert r["cos"] > 0.999
    assert r["native_pctile"] == 0.0
    for k in M.RUNGS:
        if k == "PROD":
            continue
        assert r["pref"][k] == (1.0 if orr[k] < orr["PROD"] else 0.0), k
    assert r["pref"]["circ_best"] == 1.0 and r["pref"]["RAND_SIGNED[0]"] == 0.0
    # the negated cost: every sign flips
    g, _ = _synthetic(monkeypatch, -1.0)
    monkeypatch.setattr(sys.modules[__name__], "_syn_cost", g, raising=False)
    r2 = M.meter_target("SYN", M.Cost(__name__ + ":_syn_cost"))
    for nm in M.LADDERS:
        assert abs(r2["rho"][nm] + 1.0) < 1e-12
    assert r2["cos"] < -0.999 and r2["native_pctile"] == 1.0 and r2["pref"]["circ_best"] == 0.0


def test_helpers():
    assert M.pct_in_pool(np.array([1.0, 2.0, 3.0, 3.0]), 3.0) == 0.75
    assert M.pct_in_pool(np.array([1.0, 2.0, 3.0]), 0.0) == 0.0
    assert abs(M.spearman([3, 1, 2], [3, 1, 2]) - 1.0) < 1e-12
    assert np.isnan(M.spearman([1, 1, 1], [1, 2, 3]))
    assert np.isnan(M.spearman([1, np.nan, 2], [1, 2, 3]))
    with pytest.raises(KeyError):
        M.Cost("NO_SUCH_COST")
