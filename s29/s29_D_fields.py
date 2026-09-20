#!/usr/bin/env python
"""s29/s29_D_fields.py -- lane D: the ORACLE cosine of every DISPLACEMENT FIELD the repository can
already produce and lane T's S29-L23 survey did not cover.

THE BOUND UNDER ATTACK. S29-L23: every native-free operator's output is a displacement A = c + u,
so RMSD_achievable = RMSD_prod sqrt(1 - rho_max^2) with rho the cosine between the displacement
and the direction to the native; assumption B2 is rho_max <= 0.14 (the random-shape-field
reference) for every field constructible from the present information. The bound moves if and only
if some native-free field has |cos| materially above 0.14 with the fold CI excluding it.

EVERY COSINE HERE IS ORACLE (it reads the direction to the native). The FIELDS are native-free:
each is a difference of two structures the pipeline can build without a native. Rigid-body
components are projected out of both sides (`s27/s28_A2_local.py`, S28-L23b's machinery), and the
random-shape-field reference is recomputed per target from the same code.

Fields (all native-free; production c = the DIS top-75 uniform average in its medoid frame):
  MSET_<m>     (top-m average) - c, for m in 1, 5, 10, 25, 50, 150, 250, 500 -- the m-ladder's own
               displacement, the one scalar lane T's S29-L15 says the quantum stage still chooses.
  MEDOID       (the retained set's medoid member) - c: the on-manifold end of the perception-
               distortion curve (S29-L12), which no survey rung covered.
  PROJ         (the production PROJECTION's CA trace) - c: the deployed operator that actually
               moves production 0.164 A on the reporting basis, and the only displacement in the
               shipped pipeline downstream of the average.
  CHAN_<X>     (the top-75 average under DIS+X) - c for the S27 channels X: the re-ranking
               displacement, i.e. what every alternative objective in the library buys as a field.
  CONS_TRIM    (the average of the 75 members nearest the pool medoid) - c: the consensus family's
               displacement.
  EXPAND       c scaled about its centroid to the pool's mean Rg - c: the pure de-contraction
               direction (Jensen, S29-L12), which is the field the record's contraction story
               implies and which no rung tested as a DIRECTION.

    python s29/s29_D_fields.py run [--limit N]
    python s29/s29_D_fields.py analyse

Results `s29/results/s29_D_fields_rows.jsonl`, `s29/results/s29_D_fields.json`.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
from s15 import seed as SD                 # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import s28_A2_local as A2         # noqa: E402
from s29 import s29_D_cost_audit as M      # noqa: E402

RESULTS = os.path.join(HERE, "results")
ROWS = os.path.join(RESULTS, "s29_D_fields_rows.jsonl")
M_LADDER = (1, 5, 10, 25, 50, 150, 250, 500)
CHANNELS = ("LEG", "DISTPOT", "CONTACT", "ENV", "CAGEO", "RG_LAW", "SS_MATCH", "DIS_MEAN", "CONTACT_LL")
N_REF = 16
PROD_M = 75


def field_row(pdb):
    t0 = time.time()
    from s24 import d_harness as H
    from s27 import run_pool as RP
    cand, ch, _ = RP.channels_for(pdb)
    n, k = int(cand.n), int(cand.k)
    key = RP.rng_for(pdb, "tiekey").random(k)
    dis = np.asarray(ch["DIS"], float)
    order = np.lexsort((key, dis))
    top = np.sort(order[:PROD_M])
    C0, _ = H.readout_uniform(cand, top)                                  # production
    C0 = np.asarray(C0, float)
    u = A2.remove_rigid(A2.oracle_direction(C0, cand.nat_ca), C0)         # ORACLE
    W = I.superpose_batch(np.asarray(cand.W, float), C0)                  # everything in c's frame
    row = dict(pdb=pdb, n=n, fold=int(cand.fold), fail18=bool(pdb in I.FAIL18),
               rmsd_prod=float(I.ca_rmsd(C0, cand.nat_ca)), u_rms=A2.rms(u), cos={}, disp_rms={})

    def add(name, C1):
        d = A2.remove_rigid(np.asarray(C1, float) - C0, C0)
        r = A2.rms(d)
        row["cos"][name] = A2.cosine(d, u) if r > 1e-12 else float("nan")   # ORACLE
        row["disp_rms"][name] = float(r)
    for m in M_LADDER:
        Cm, _ = H.readout_uniform(cand, np.sort(order[:m]))
        add(f"MSET_{m}", Cm)
    P = I.pairwise_rmsd(W[top])
    add("MEDOID", W[top][int(I.medoid(P))])
    pr = I.project(C0, cand.seq, cand.fold)
    add("PROJ", np.asarray(pr["ca"], float))
    for X in CHANNELS:
        if X not in ch:
            continue
        E = RP.combine(ch, ["DIS", X])
        oX = np.lexsort((key, np.asarray(E, float)))
        CX, _ = H.readout_uniform(cand, np.sort(oX[:PROD_M]))
        add(f"CHAN_{X}", CX)
    # the consensus trim: the 75 members nearest the pool medoid
    Pall = I.pairwise_rmsd(W)
    cons = Pall.sum(1) / max(k - 1, 1)
    oc = np.lexsort((key, cons))
    Cc, _ = H.readout_uniform(cand, np.sort(oc[:PROD_M]))
    add("CONS_TRIM", Cc)
    # pure de-contraction: scale c about its centroid to the pool's mean Rg
    rg_pool = float(np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1)).mean())
    cen = C0.mean(0)
    rg0 = float(np.sqrt(((C0 - cen) ** 2).sum(1).mean()))
    add("EXPAND", cen + (C0 - cen) * (rg_pool / max(rg0, 1e-12)))
    rng = SD.stable_rng(pdb, "s28A2_cosref", salt=A2.SALT)                 # A2's own reference draws
    row["cos_random_ref"] = [A2.cosine(A2.remove_rigid(rng.normal(size=C0.shape), C0), u) for _ in range(N_REF)]
    row["secs"] = time.time() - t0
    return row


def run(pdbs):
    os.makedirs(RESULTS, exist_ok=True)
    done = set()
    if os.path.exists(ROWS):
        with open(ROWS, encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    done.add(json.loads(line)["pdb"])
    t0 = time.time()
    for q, pdb in enumerate(pdbs):
        if pdb in done:
            continue
        row = field_row(pdb)
        with open(ROWS, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(row) + "\n")
        best = max(row["cos"], key=lambda kk: abs(row["cos"][kk]) if np.isfinite(row["cos"][kk]) else -1)
        print(f"  [{q+1}/{len(pdbs)}] {pdb} best |cos| {best} {row['cos'][best]:+.3f} | PROJ {row['cos'].get('PROJ', float('nan')):+.3f} "
              f"EXPAND {row['cos'].get('EXPAND', float('nan')):+.3f} MSET_1 {row['cos'].get('MSET_1', float('nan')):+.3f} "
              f"{row['secs']:.1f}s (elapsed {(time.time()-t0)/60:.1f} min)", flush=True)
    print("done:", ROWS)


def second_order(pdbs):
    """S29-L23 assumption B3: the bound is stated first order in the step size, with the rigid
    body projected out.  For each field, take its OWN best step (the ORACLE s* that minimises the
    RMSD along it) and compare the bound's predicted RMSD, RMSD_prod * sqrt(1 - rho^2), with the
    RMSD actually measured after Kabsch re-alignment.  The gap is the neglected term AT THE
    LARGEST STEP A REAL ARM TAKES, not at s = 0.14."""
    rows = []
    for q, pdb in enumerate(pdbs):
        from s24 import d_harness as H
        from s27 import run_pool as RP
        cand, ch, _ = RP.channels_for(pdb)
        n, k = int(cand.n), int(cand.k)
        key = RP.rng_for(pdb, "tiekey").random(k)
        dis = np.asarray(ch["DIS"], float)
        order = np.lexsort((key, dis))
        top = np.sort(order[:PROD_M])
        C0, _ = H.readout_uniform(cand, top); C0 = np.asarray(C0, float)
        u = A2.remove_rigid(A2.oracle_direction(C0, cand.nat_ca), C0)          # ORACLE
        W = I.superpose_batch(np.asarray(cand.W, float), C0)
        rec = dict(pdb=pdb, n=n, rmsd_prod=float(I.ca_rmsd(C0, cand.nat_ca)), fields={})
        fields = {}
        for m in (25, 150, 500):
            Cm, _ = H.readout_uniform(cand, np.sort(order[:m]))
            fields[f"MSET_{m}"] = np.asarray(Cm, float)
        pr = I.project(C0, cand.seq, cand.fold)
        fields["PROJ"] = np.asarray(pr["ca"], float)
        cen = C0.mean(0); rg0 = float(np.sqrt(((C0 - cen) ** 2).sum(1).mean()))
        rgp = float(np.sqrt(((W - W.mean(1, keepdims=True)) ** 2).sum(2).mean(1)).mean())
        fields["EXPAND"] = cen + (C0 - cen) * (rgp / max(rg0, 1e-12))
        for nm, C1 in fields.items():
            d = A2.remove_rigid(np.asarray(C1, float) - C0, C0)
            if A2.rms(d) < 1e-12:
                continue
            rho = A2.cosine(d, u)
            dh = d / A2.rms(d)
            #: the ORACLE best step along this field, on a fine grid (the largest step a real arm
            #: could take, not a nominal one), and the RMSD actually measured there
            ss = np.linspace(-3.0, 3.0, 241)
            r = np.array([I.ca_rmsd(C0 + t * dh, cand.nat_ca) for t in ss])     # ORACLE
            b = int(np.argmin(r))
            pred = rec["rmsd_prod"] * float(np.sqrt(max(1 - rho ** 2, 0.0)))
            rec["fields"][nm] = dict(cos=rho, disp_rms=A2.rms(d), best_step=float(ss[b]),
                                     rmsd_at_best=float(r[b]), predicted_by_bound=pred,
                                     residual=float(r[b] - pred),
                                     rel_residual=float((r[b] - pred) / max(pred, 1e-12)))
        rows.append(rec)
        if (q + 1) % 6 == 0:
            print(f"  [{q+1}/{len(pdbs)}]", flush=True)
    out = dict(check="S29-L23 assumption B3: the neglected term at the ORACLE best step", n=len(rows), rows=rows)
    names = sorted({k for r in rows for k in r["fields"]})
    print(f"  {'field':12s} {'mean cos':>9s} {'best step':>10s} {'measured':>9s} {'bound pred':>11s} {'residual':>9s} {'rel':>7s}")
    for nm in names:
        v = [r["fields"][nm] for r in rows if nm in r["fields"]]
        if not v:
            continue
        out[nm] = dict(n=len(v), mean_cos=float(np.mean([x["cos"] for x in v])),
                       mean_best_step=float(np.mean([x["best_step"] for x in v])),
                       mean_rmsd_at_best=float(np.mean([x["rmsd_at_best"] for x in v])),
                       mean_predicted=float(np.mean([x["predicted_by_bound"] for x in v])),
                       mean_residual=float(np.mean([x["residual"] for x in v])),
                       mean_rel_residual=float(np.mean([x["rel_residual"] for x in v])))
        d = out[nm]
        print(f"  {nm:12s} {d['mean_cos']:+9.4f} {d['mean_best_step']:+10.3f} {d['mean_rmsd_at_best']:9.4f} "
              f"{d['mean_predicted']:11.4f} {d['mean_residual']:+9.4f} {100*d['mean_rel_residual']:+6.1f}%")
    ST.save_atomic(os.path.join(RESULTS, "s29_D_fields_b3.json"), out, module_file=__file__)
    return out


def analyse():
    rows = [json.loads(l) for l in open(ROWS, encoding="utf-8") if l.strip()]
    rows.sort(key=lambda r: r["pdb"])
    pdbs = [r["pdb"] for r in rows]
    folds = ST.pinned_folds(pdbs)
    fail = np.array([r["fail18"] for r in rows])
    names = sorted({k for r in rows for k in r["cos"]})
    ref = float(np.mean([np.mean(np.abs(r["cos_random_ref"])) for r in rows]))
    out = dict(check="ORACLE cosines of native-free displacement fields (attack on S29-L23 assumption B2)",
               n=len(rows), random_ref_mean_abs=ref, prod_rmsd=float(np.mean([r["rmsd_prod"] for r in rows])),
               fields={})
    print(f"  random-shape-field reference |cos| = {ref:.3f} (the bound's B2 value is 0.14); n = {len(rows)}")
    print(f"  {'field':14s} {'mean cos':>9s} {'fold CI':>20s} {'|cos| mean':>10s} {'disp RMS':>9s} {'FAIL18':>8s} {'108':>8s} "
          f"{'implied RMSD':>12s}  verdict")
    for nm in names:
        c = np.array([r["cos"].get(nm, np.nan) for r in rows], float)
        ok = np.isfinite(c)
        if ok.sum() < 10:
            continue
        cmp0 = ST.compare(c[ok], np.zeros(int(ok.sum())), folds[ok], label=f"cos({nm}, u) vs 0", seed_parts=("s29Dfield",))
        absm = float(np.abs(c[ok]).mean())
        rmsd_p = float(np.mean([r["rmsd_prod"] for r in rows]))
        implied = rmsd_p * float(np.sqrt(max(1.0 - cmp0["effect"] ** 2, 0.0)))
        beats = bool(abs(cmp0["effect"]) - 1.96 * cmp0["se"] > ref)
        out["fields"][nm] = dict(mean=cmp0["effect"], se=cmp0["se"], ci95_fold=cmp0["ci95_fold"],
                                 ci95_iid=cmp0["ci95_iid"], median=cmp0["median_effect"], n=int(ok.sum()),
                                 mean_abs=absm, disp_rms=float(np.mean([r["disp_rms"].get(nm, np.nan) for r in rows])),
                                 fail18=float(np.nanmean(c[fail & ok])), other=float(np.nanmean(c[(~fail) & ok])),
                                 implied_rmsd_at_best_step=implied, beats_random_reference=beats,
                                 n_positive=int((c[ok] > 0).sum()))
        v = out["fields"][nm]
        print(f"  {nm:14s} {v['mean']:+9.4f} [{v['ci95_fold'][0]:+.3f},{v['ci95_fold'][1]:+.3f}] {absm:10.3f} "
              f"{v['disp_rms']:9.3f} {v['fail18']:+8.3f} {v['other']:+8.3f} {implied:12.4f}  "
              f"{'BEATS the reference' if beats else 'inside the reference'}")
    best = max(out["fields"], key=lambda kk: abs(out["fields"][kk]["mean"]))
    out["best_field"] = best
    out["best_mean_cos"] = out["fields"][best]["mean"]
    out["any_beats_reference"] = [kk for kk, v in out["fields"].items() if v["beats_random_reference"]]
    out["verdict"] = ("B2 SURVIVES this survey: no field's mean |cos| clears the random reference %.3f "
                      "with a 2-SE margin (best %s at %+.4f)" % (ref, best, out["fields"][best]["mean"])) \
        if not out["any_beats_reference"] else \
        ("B2 MOVES: %s clear the random reference" % ", ".join(out["any_beats_reference"]))
    print(" ", out["verdict"])
    ST.save_atomic(os.path.join(RESULTS, "s29_D_fields.json"), out, module_file=__file__)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["run", "analyse", "b3"])
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    from s25 import phys_lib as P
    pdbs = P.targets()[:a.limit] if a.limit else P.targets()
    if a.mode == "run":
        run(pdbs)
    elif a.mode == "b3":
        second_order(pdbs)
    else:
        analyse()


if __name__ == "__main__":
    main()
