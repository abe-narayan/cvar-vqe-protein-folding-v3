"""s16/verify_stats.py -- VERIFY workstream. Re-analysis of the two artefact-only claims.

Nothing here re-runs science. It reads `s16/results/integrate.json` and
`s16/results/diversity.json` and re-does the statistics with the TARGET as the unit of
analysis, which the BRIEF (3.4) makes binding and which `integrate._boot` /
`diversity._boot` do not do -- both resample ROWS, and a row is a (target x seed x
generator) cell.

Outputs `s16/results/verify_stats.json`.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
RESULTS = os.path.join(HERE, "results")

from s15 import seed as SD  # noqa: E402


def boot_rows(d, rng, B=20000):
    """the coordinator's convention: resample ROWS i.i.d."""
    d = np.asarray(d, float)
    k = len(d)
    m = d[rng.integers(0, k, size=(B, k))].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)), k


def boot_targets(pairs, rng, B=20000):
    """the BRIEF's convention: the TARGET is the unit.

    `pairs` is [(pdb, value)].  Per-target means first, then a cluster bootstrap over
    targets (resample targets with replacement, take the mean of their per-target means).
    """
    byp = {}
    for p, v in pairs:
        byp.setdefault(p, []).append(float(v))
    ps = sorted(byp)
    mu = np.array([np.mean(byp[p]) for p in ps])
    k = len(ps)
    bm = mu[rng.integers(0, k, size=(B, k))].mean(1)
    w = int((mu > 0).sum())
    l = int((mu < 0).sum())
    return (float(mu.mean()), float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5)),
            k, float(np.median(mu)), w, l, {p: float(np.mean(byp[p])) for p in ps})


def integrate_reanalysis(rng):
    rows = json.load(open(os.path.join(RESULTS, "integrate.json")))["rows"]
    MS = sorted({r["m"] for r in rows})
    gens = ["cvar0.1", "cvar0.25", "cvar0.5", "cvar1.0", "bestofN", "uniform"]
    out = {"n_rows": len(rows), "n_targets": len({r["pdb"] for r in rows}), "effects": {}}

    contrasts = {
        "LEGACY = leg - ctrl": lambda r: r["leg"] - r["ctrl"],
        "AMBER  = amb - ctrl": lambda r: r["amb"] - r["ctrl"],
        "AMBER  = amb - rnd": lambda r: r["amb"] - r["rnd"],
        "LEGACY = leg - rnd": lambda r: r["leg"] - r["rnd"],
        "INTER": lambda r: (r["legamb"] - r["amb"]) - (r["leg"] - r["ctrl"]),
        "tors - ctrl": lambda r: r["tors"] - r["ctrl"],
    }
    print("=" * 96)
    print("INTEGRATE -- rows-as-unit (coordinator) vs TARGET-as-unit (BRIEF 3.4)")
    print("=" * 96)
    for m in MS:
        sel = [r for r in rows if r["m"] == m]
        print(f"\n--- m = {m}   ({len(sel)} rows, {len({r['pdb'] for r in sel})} targets)")
        print(f"  {'contrast':<22}{'ROWS mean [95% CI]':>30}{'TARGET mean [95% CI]':>30}"
              f"{'med':>8}{'W/L':>8}{'CI x':>7}")
        for cname, fn in contrasts.items():
            a, lo, hi, nr = boot_rows([fn(r) for r in sel], rng)
            ta, tlo, thi, nt, tmed, w, l, per = boot_targets([(r["pdb"], fn(r)) for r in sel], rng)
            widen = (thi - tlo) / max(hi - lo, 1e-12)
            print(f"  {cname:<22}{a:>+11.4f} [{lo:+.4f},{hi:+.4f}]"
                  f"{ta:>+11.4f} [{tlo:+.4f},{thi:+.4f}]{tmed:>+8.3f}{f'{w}/{l}':>8}{widen:>7.2f}")
            out["effects"].setdefault(str(m), {})[cname] = {
                "rows": {"mean": a, "ci95": [lo, hi], "n": nr},
                "target": {"mean": ta, "ci95": [tlo, thi], "n": nt, "median": tmed,
                           "wl": [w, l], "per_target": per},
                "ci_widen": widen}

    # ---- is the m=75 AMBER effect carried by one or two targets? leave-one-target-out
    print("\n" + "=" * 96)
    print("AMBER (amb - ctrl) at m = 75: leave-one-TARGET-out, and per generator")
    print("=" * 96)
    sel = [r for r in rows if r["m"] == 75]
    ta, tlo, thi, nt, tmed, w, l, per = boot_targets(
        [(r["pdb"], r["amb"] - r["ctrl"]) for r in sel], rng)
    print("  per target: " + "  ".join(f"{p} {v:+.3f}" for p, v in sorted(per.items())))
    lo1 = {}
    for p in sorted(per):
        s2 = [r for r in sel if r["pdb"] != p]
        a2, l2, h2, *_ = boot_targets([(r["pdb"], r["amb"] - r["ctrl"]) for r in s2], rng)[:3]
        lo1[p] = [a2, l2, h2]
        print(f"  drop {p}: {a2:+.4f} [{l2:+.4f},{h2:+.4f}]"
              + ("   <-- CI now contains zero" if l2 <= 0 <= h2 else ""))
    out["amber75_loo"] = lo1
    out["amber75_per_target"] = per

    print("\n  per generator (target as unit, n = 9):")
    pg = {}
    for g in gens:
        s2 = [r for r in sel if r["gen"] == g]
        a2, l2, h2, n2, md2, w2, l2c, _ = boot_targets(
            [(r["pdb"], r["amb"] - r["ctrl"]) for r in s2], rng)
        pg[g] = {"mean": a2, "ci95": [l2, h2], "median": md2, "wl": [w2, l2c]}
        print(f"    {g:<10}{a2:>+9.4f} [{l2:+.4f},{h2:+.4f}]  median {md2:+.3f}  {w2}W/{l2c}L")
    out["amber75_per_gen"] = pg

    # ---- amb vs rnd at m=75, the comparison the FINDINGS quote
    print("\n  amb - rnd at m = 75 (the quoted contrast), per generator, target as unit:")
    pg2 = {}
    for g in gens:
        s2 = [r for r in sel if r["gen"] == g]
        a2, l2, h2, n2, md2, w2, l2c, _ = boot_targets(
            [(r["pdb"], r["amb"] - r["rnd"]) for r in s2], rng)
        pg2[g] = {"mean": a2, "ci95": [l2, h2], "median": md2, "wl": [w2, l2c]}
        print(f"    {g:<10}{a2:>+9.4f} [{l2:+.4f},{h2:+.4f}]  median {md2:+.3f}  {w2}W/{l2c}L")
    out["amber75_vs_rnd_per_gen"] = pg2

    # ---- seed-level replication: does the m=75 AMBER effect hold in each seed alone?
    print("\n  by SEED (target as unit within each seed), amb - ctrl at m = 75:")
    bs = {}
    for s in sorted({r["seed"] for r in sel}):
        s2 = [r for r in sel if r["seed"] == s]
        a2, l2, h2, n2, md2, w2, l2c, _ = boot_targets(
            [(r["pdb"], r["amb"] - r["ctrl"]) for r in s2], rng)
        bs[s] = {"mean": a2, "ci95": [l2, h2]}
        print(f"    seed {s}: {a2:>+9.4f} [{l2:+.4f},{h2:+.4f}]  {w2}W/{l2c}L")
    out["amber75_by_seed"] = bs

    # ---- INTEGRATE_FINDINGS section 3: every CVaR-VQE arm loses to uniform / bestofN
    print("\n" + "=" * 96)
    print("INTEGRATE_FINDINGS section 3, re-read with the TARGET as the unit (n = 9)")
    print("=" * 96)
    vq = {}
    for m in MS:
        print(f"\n  --- m = {m}   (positive = the VQE arm is WORSE)")
        print(f"  {'arm':<12}{'minus bestofN (rows)':>28}{'minus bestofN (TARGET)':>30}"
              f"{'W/L':>7}")
        for g in ("cvar0.1", "cvar0.25", "cvar0.5", "cvar1.0"):
            for ctl in ("bestofN", "uniform"):
                pr = []
                for p in sorted({r["pdb"] for r in rows}):
                    for s in sorted({r["seed"] for r in rows}):
                        a = [r for r in rows if r["pdb"] == p and r["seed"] == s
                             and r["m"] == m and r["gen"] == g]
                        b = [r for r in rows if r["pdb"] == p and r["seed"] == s
                             and r["m"] == m and r["gen"] == ctl]
                        if a and b:
                            pr.append((p, a[0]["ctrl"] - b[0]["ctrl"]))
                ra, rlo, rhi, _ = boot_rows([v for _p, v in pr], rng)
                ta, tlo, thi, nt, tmed, w, l, _ = boot_targets(pr, rng)
                vq[f"{g}-{ctl}-m{m}"] = {
                    "rows": [ra, rlo, rhi], "target": [ta, tlo, thi],
                    "median": tmed, "wl": [w, l]}
                if ctl == "bestofN":
                    line = (f"  {g:<12}{ra:>+11.3f} [{rlo:+.3f},{rhi:+.3f}]"
                            f"{ta:>+13.3f} [{tlo:+.3f},{thi:+.3f}]{f'{w}/{l}':>7}")
                else:
                    line = (f"  {'  vs unif':<12}{ra:>+11.3f} [{rlo:+.3f},{rhi:+.3f}]"
                            f"{ta:>+13.3f} [{tlo:+.3f},{thi:+.3f}]{f'{w}/{l}':>7}")
                print(line)
    out["vqe_controls"] = vq
    return out


def diversity_reanalysis(rng):
    rows = json.load(open(os.path.join(RESULTS, "diversity.json")))["rows"]
    gens = ["cvar0.1", "cvar0.25", "cvar0.5", "cvar1.0", "bestofN", "uniform"]
    out = {"n_rows": len(rows)}

    print("\n" + "=" * 96)
    print("DIVERSITY -- is the Krogh-Vedelsby residual a measurement or a theorem check?")
    print("=" * 96)
    res = np.array([r["readout"] - r["kv_pred"] for r in rows])
    # the identity, in the squared form that is the actual algebra
    id_res = np.array([r["readout"] ** 2 - (r["err2_common_frame"] - r["div2"]) for r in rows])
    print(f"  residual (RMSD units)     mean {res.mean():+.3e}  max|.| {np.abs(res).max():.3e}")
    print(f"  residual (squared form)   mean {id_res.mean():+.3e}  max|.| {np.abs(id_res).max():.3e}")
    scale = np.array([r["readout"] ** 2 for r in rows])
    print(f"  relative to readout^2     max|rel| {np.abs(id_res / scale).max():.3e}"
          "   <-- machine epsilon is ~2.2e-16")
    out["kv_residual"] = {"rmsd_units_max_abs": float(np.abs(res).max()),
                          "squared_form_max_abs": float(np.abs(id_res).max()),
                          "max_relative": float(np.abs(id_res / scale).max())}

    print("\n  THE TERM THE IDENTITY ACTUALLY USES is `err2_common_frame`, NOT the free-")
    print("  superposition `mean_member_quad` that the findings table prints. Both, at m=75:")
    print(f"  {'generator':<12}{'readout':>9}{'free-sup RMSq':>15}{'COMMON-FRAME err':>18}"
          f"{'diversity':>11}")
    tab = {}
    for g in gens:
        sel = [r for r in rows if r["gen"] == g and r["m"] == 75]
        ro = np.mean([r["readout"] for r in sel])
        fq = np.mean([r["mean_member_quad"] for r in sel])
        cf = np.mean([np.sqrt(r["err2_common_frame"]) for r in sel])
        dv = np.mean([r["diversity"] for r in sel])
        tab[g] = {"readout": ro, "free_sup_quad": fq, "common_frame_err": cf, "diversity": dv}
        print(f"  {g:<12}{ro:>9.3f}{fq:>15.3f}{cf:>18.3f}{dv:>11.3f}")
    sp = lambda k: (max(v[k] for v in tab.values()) - min(v[k] for v in tab.values()))  # noqa
    print(f"  spread:      readout {sp('readout'):.3f}   free-sup {sp('free_sup_quad'):.3f}   "
          f"COMMON-FRAME {sp('common_frame_err'):.3f}   diversity {sp('diversity'):.3f}")
    out["m75_table"] = tab

    # ---- attribution: decompose the readout^2 spread into error and diversity parts,
    #      PER TARGET, against the `uniform` generator as the reference.
    print("\n  ATTRIBUTION, per target, m = 75, each generator vs `uniform`:")
    print("     d(readout^2) = d(err2_common_frame) - d(div2)   [exact]")
    print(f"  {'generator':<12}{'d readout^2':>13}{'from ERROR':>13}{'from DIVERSITY':>16}"
          f"{'error share':>13}{'targets':>9}")
    att = {}
    for g in gens:
        if g == "uniform":
            continue
        de, dd, dr = [], [], []
        for p in sorted({r["pdb"] for r in rows}):
            A = [r for r in rows if r["pdb"] == p and r["m"] == 75 and r["gen"] == g]
            B = [r for r in rows if r["pdb"] == p and r["m"] == 75 and r["gen"] == "uniform"]
            de.append(np.mean([r["err2_common_frame"] for r in A]) -
                      np.mean([r["err2_common_frame"] for r in B]))
            dd.append(-(np.mean([r["div2"] for r in A]) - np.mean([r["div2"] for r in B])))
            dr.append(np.mean([r["readout"] ** 2 for r in A]) -
                      np.mean([r["readout"] ** 2 for r in B]))
        de = np.array(de); dd = np.array(dd); dr = np.array(dr)
        share = abs(de.mean()) / max(abs(de.mean()) + abs(dd.mean()), 1e-12)
        att[g] = {"d_readout2": float(dr.mean()), "from_error": float(de.mean()),
                  "from_diversity": float(dd.mean()), "error_share": float(share),
                  "per_target_error": de.tolist(), "per_target_div": dd.tolist()}
        print(f"  {g:<12}{dr.mean():>+13.3f}{de.mean():>+13.3f}{dd.mean():>+16.3f}"
              f"{share:>13.2f}{len(de):>9}")
    out["attribution_vs_uniform"] = att

    # ---- the claim "near-equal mean member error", tested with the target as the unit
    print("\n  IS THE MEMBER ERROR EQUAL ACROSS GENERATORS?  paired per target, m = 75,")
    print("  each generator minus `uniform`, target as the unit (n = 9):")
    print(f"  {'generator':<12}{'free-sup RMSq':>26}{'COMMON-FRAME err':>28}")
    eq = {}
    for g in gens:
        if g == "uniform":
            continue
        pf, pc = [], []
        for p in sorted({r["pdb"] for r in rows}):
            A = [r for r in rows if r["pdb"] == p and r["m"] == 75 and r["gen"] == g]
            B = [r for r in rows if r["pdb"] == p and r["m"] == 75 and r["gen"] == "uniform"]
            pf.append((p, np.mean([r["mean_member_quad"] for r in A]) -
                       np.mean([r["mean_member_quad"] for r in B])))
            pc.append((p, np.mean([np.sqrt(r["err2_common_frame"]) for r in A]) -
                       np.mean([np.sqrt(r["err2_common_frame"]) for r in B])))
        a1, l1, h1, *_ = boot_targets(pf, rng)[:3]
        a2, l2, h2, *_ = boot_targets(pc, rng)[:3]
        eq[g] = {"free_sup": [a1, l1, h1], "common_frame": [a2, l2, h2]}
        print(f"  {g:<12}{a1:>+10.3f} [{l1:+.3f},{h1:+.3f}]{a2:>+12.3f} [{l2:+.3f},{h2:+.3f}]")
    out["member_error_equality"] = eq
    return out


def main():
    rng = SD.stable_rng("s16verify", "stats")
    out = {"integrate": integrate_reanalysis(rng), "diversity": diversity_reanalysis(rng)}
    with open(os.path.join(RESULTS, "verify_stats.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote s16/results/verify_stats.json")


if __name__ == "__main__":
    main()
