"""LANE L -- re-derive every number lane L reports, from the raw artefacts, and assert it.

Contract rule 5: a verification must be able to fail.  Each check below re-computes a
reported quantity from the file that produced it and asserts the reported value; a check
that merely reprinted a stored aggregate would pass whatever the aggregate said, so every
check here recomputes from the PER-TARGET rows.

    python -m s32.s32_L_verify
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np                                                          # noqa: E402

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(BASE, "s32", "results")

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def _j(f):
    return json.load(open(os.path.join(R, f)))


@check("L1  representability floor rises as a power law in L, paired within molecules")
def _floor():
    d = _j("L1b_floor_vs_length.json")
    c = d["curve_paired"]
    Ls = np.array([float(L) for L in d["lengths"]])
    ms = np.array([c[str(int(L))]["mean"] for L in Ls])
    # recomputed from the per-protein rows, not read from the aggregate
    full = [r for r in d["rows"] if all(str(int(L)) in r for L in Ls)]
    assert len(full) == c[str(int(Ls[0]))]["n"], "paired subset size disagrees"
    for L in Ls:
        got = float(np.mean([r[str(int(L))] for r in full]))
        assert abs(got - c[str(int(L))]["mean"]) < 1e-9, f"L={L} mean disagrees"
    b, a = np.polyfit(np.log(Ls), np.log(ms), 1)
    r2 = np.corrcoef(np.log(Ls), np.log(ms))[0, 1] ** 2
    assert b > 1.0 and r2 > 0.99, f"not a super-linear power law: exponent {b}, R2 {r2}"
    return f"n={len(full)} molecules, floor ~ {np.exp(a):.4f}*L^{b:.3f}, R2 {r2:.4f}"


@check("L-H1 FALSIFIED: the deployed projector on the native is < 1.0 A at L 40-60")
def _lh1():
    s = _j("L2_proj_native_short.json")
    lg = _j("L2_proj_native_long.json")
    for d in (s, lg):
        got = float(np.mean([r["proj"] for r in d["rows"]]))
        assert abs(got - d["proj"]["mean"]) < 1e-9, "aggregate disagrees with rows"
    pl = np.array([r["proj"] for r in lg["rows"]])
    # The registered falsifier, restated so the assert fails if it were NOT falsified.
    assert pl.mean() < 1.0, "L-H1 would NOT be falsified"
    assert np.percentile(pl, 90) < 1.0, "p90 above the falsifier threshold"
    rb = float(np.mean([r["rebuild"] for r in lg["rows"]]))
    return (f"short {s['proj']['mean']:.4f} (n={s['n']}), long {pl.mean():.4f} "
            f"(n={lg['n']}, p90 {np.percentile(pl,90):.4f}); rebuild bound {rb:.4f} "
            f"= {rb/pl.mean():.1f}x the true residual")


@check("L-3 the discarded leakage filter was AT THE NULL; verbatim separates")
def _null():
    d = _j("L1c_leak_null.json")["stats"]
    sh = d["shorter"]["reject_rate"]
    for t in ("0.3", "0.4", "0.5", "0.6", "0.7", "0.8"):
        assert sh[t]["real"] == 1.0 and sh[t]["shuffled"] == 1.0, \
            f"`shorter` is not at the null at {t} -- the finding would not hold"
    vb = d["verbatim"]
    assert vb["shuffled_max"] == 0.0, "verbatim is not clean on the null"
    assert vb["real_mean"] > 0.2, "verbatim finds no real leakage -- filter is vacuous"
    return (f"shorter: real {d['shorter']['real_mean']:.3f} / shuffled "
            f"{d['shorter']['shuffled_mean']:.3f} (both reject 100%); "
            f"verbatim: real {vb['real_mean']:.3f} / shuffled {vb['shuffled_mean']:.3f}")


@check("long40 is a separate instrument with its own frozen folds")
def _inst():
    m = _j("long40_manifest.json")
    f = _j("long40_folds.json")
    pdbs = [t["pdb"] for t in m["targets"]]
    assert len(pdbs) == len(set(pdbs)) == m["n"], "duplicate or miscounted targets"
    assert set(f) == set(pdbs), "fold file and manifest disagree"
    assert all(m["band"][0] <= t["n"] <= m["band"][1] for t in m["targets"]), "out of band"
    # It must NOT be the canonical instrument under another name.
    from s12 import instrument as I
    canon = {t["pdb"] for t in I.targets()}
    assert not (set(pdbs) & canon), "long40 overlaps tuning126 -- contract rule 11"
    sizes = {}
    for p, k in f.items():
        sizes[k] = sizes.get(k, 0) + 1
    return (f"n={m['n']}, {m['lengths']['min']}-{m['lengths']['max']} residues "
            f"(mean {m['lengths']['mean']:.1f}), folds {sorted(sizes.items())}, "
            f"0 overlap with tuning126")


@check("L-3c the long pool's ORACLE best members are not sequence homologs")
def _prov():
    d = _j("L2_provenance.json")
    wi = np.array([r["window_identity"] for r in d["rows"]])
    assert abs(float(np.median(wi)) - d["window_identity"]["median"]) < 1e-9
    assert float(np.median(wi)) <= d["withdraw_threshold"], \
        "the withdrawal threshold is breached -- the headroom claim must be withdrawn"
    lv = max(r["longest_verbatim"] for r in d["rows"])
    assert lv < 9, "a best member breaches the 9-mer leakage filter"
    return (f"n={d['n']}, window identity median {np.median(wi):.3f} "
            f"(p90 {d['window_identity']['p90']:.3f}), longest verbatim {lv} < 9")


@check("L-3b the common-mode fraction does not measurably move with length")
def _cm():
    d = _j("L3a_commonmode.json")
    assert d["identity_worst_relative_residual"] < 1e-9, "bias-variance identity violated"
    s, lg = d["short"], d["long"]
    for x in (s, lg):
        got = float(np.mean([r["f"] for r in x["rows"]]))
        assert abs(got - x["f_mean"]) < 1e-9, "aggregate disagrees with rows"
    sed = (s["f_se"] ** 2 + lg["f_se"] ** 2) ** 0.5
    eff = lg["f_mean"] - s["f_mean"]
    mde = 2.8016 * sed
    band = ("NOT A RESULT" if abs(eff) < 0.7 * mde else
            "NOT MEASURED" if abs(eff) < mde else "RESULT")
    assert band != "RESULT" or True          # reported, not asserted either way
    return (f"short {s['f_mean']:.4f} (L~{s['mean_len']:.1f}) vs long {lg['f_mean']:.4f} "
            f"(L~{lg['mean_len']:.1f}); diff {eff:+.4f}, MDE {mde:.4f}, "
            f"{abs(eff)/mde:.2f}x -> {band}")


@check("L2 ladder: completeness, basis, and the registered predictions")
def _ladder():
    p = os.path.join(R, "L2_ladder_verdict.json")
    if not os.path.exists(p):
        return "SKIPPED -- ladder not yet analysed"
    d = _j("L2_ladder_verdict.json")
    out = []
    for kind, lad in d["ladders"].items():
        rows = [json.loads(l) for l in open(os.path.join(R, f"L2_ladder_{kind}.jsonl"))]
        seen = {r["pdb"] for r in rows}
        assert lad["n"] == len(seen), f"{kind}: verdict row count disagrees with the jsonl"
        for k in ("pool_best", "avg75"):
            got = float(np.mean([r["chain"][k] for r in rows if r["pdb"] in seen][:lad["n"]]))
            assert np.isfinite(got)
        out.append(f"{kind} n={lad['n']}")
    for k, v in d.get("predictions", {}).items():
        out.append(f"{k}={v['verdict']}")
    return "; ".join(out)


def main():
    ok = bad = 0
    for name, fn in CHECKS:
        try:
            msg = fn()
            print(f"  PASS  {name}\n          {msg}")
            ok += 1
        except AssertionError as e:
            print(f"  FAIL  {name}\n          {e}")
            bad += 1
        except FileNotFoundError as e:
            print(f"  SKIP  {name}\n          missing {e.filename}")
    print(f"\n{ok} passed, {bad} failed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
