"""s16/repair_report.py -- the REPAIR workstream's analysis.  No AMBER runs here.

Reads `s16/results/repair_sweep.json` (written by `s16/repair.py`) and produces:

  Q1  the (k, depth) table -- ACCURACY and VALIDITY as SEPARATE axes, never substituted for
      one another -- with the convergence gate applied, the excluded count printed at every
      setting, the ROTATED-LAB-FRAME NULL re-measured at every setting (it is zero by rigid
      invariance, so a setting whose null is non-zero is NOT REPORTABLE), and the effect on
      the frame-reproducible stratum stated separately.

  Q2  the GEOMETRY of what the two repair operators actually do to the structure: the cosine
      of their displacement with the true residual (ORACLE), its magnitude relative to the
      residual, and where it sits in the superposed-CA Jacobian spectrum -- with the
      ZERO-INFORMATION (constant alpha-helix) and RANDOM-DIRECTION-at-matched-magnitude
      controls that AUDIT showed this programme's evidence needs.

Statistics: TARGET is the unit, n = 126.  `I.paired`'s i.i.d. paired bootstrap is the
headline, with a FOLD-CLUSTERED twin printed beside it (5 clusters, so it is the weaker
interval, not the stronger one -- AUDIT's convention, kept).  Median and W/L beside every
mean.  The AMBER path has NO RNG and reproduces bit-identically, so the 0.08 A multi-start
false-positive floor DOES NOT apply; the floor that does is this path's own, ~0.004 A gated
and up to 0.038 A ungated.

    python -m s16.repair_report
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import json
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
RESULTS = os.path.join(HERE, "results")

from s12 import instrument as I                            # noqa: E402
from s15 import align_lib as AL                            # noqa: E402
from s15.seed import stable_rng                            # noqa: E402
from s14.avgspace import top75_windows                     # noqa: E402
from s14.ener_avgrefine import _kabsch, torsions_of        # noqa: E402
from s15.phys_repl import averaged_backbone_from, conc_verdict   # noqa: E402
from s16.repair import (SETTINGS, BY_NAME, PASSES, STABLE_DE, N_RAND, _sup,
                        shape_subsample)                   # noqa: E402
from core.amber import CONVERGE_MAX_KCAL                   # noqa: E402
import peptide_db as PDB                                   # noqa: E402

#: the standing frame-null PASS band, taken unchanged from `s16/energy_lib.py` so the two
#: workstreams' nulls are judged by one rule.
FRAME_TOL_MEAN, FRAME_TOL_MAX = 0.005, 0.05


def load(passes=("A", "B", "C")):
    """Merge the passes.  Pass A carries the base block (input CA, native CA, the raw
    average, the zero-information helix, the ideal-geometry projection); later passes carry
    only their own settings and are merged onto it BY PDB.  A target is kept only if every
    requested pass produced it, and the count is printed."""
    base = None
    have = []
    for p in passes:
        fp = os.path.join(RESULTS, f"repair_{p}.json")
        if not os.path.exists(fp):
            continue
        with open(fp) as fh:
            d = json.load(fh)
        have.append((p, {r["pdb"]: r for r in d["per_target"]}, d))
    if not have or have[0][0] != "A":
        raise SystemExit("pass A (the base) is required")
    base = have[0][1]
    keep = set(base)
    for p, m, _ in have[1:]:
        keep &= set(m)
    rows = []
    for pdbid in sorted(keep):
        r = dict(base[pdbid])
        for p, m, _ in have[1:]:
            for k, v in m[pdbid].items():
                if k.startswith("f0_") or k.startswith("f1_"):
                    r[k] = v
        rows.append(r)
    names = [nm for p, _, d in have for nm, _, _ in [tuple(x) for x in d["settings"]]]
    print("passes merged: " + ", ".join(f"{p}(n={len(m)})" for p, m, _ in have)
          + f"  -> common n = {len(rows)};  settings: {sorted(set(names))}")
    return sorted(set(names)), rows


# ------------------------------------------------------------------ statistics
def fold_boot(d, folds, n_boot=4000, seed=0):
    """Fold-CLUSTERED paired bootstrap.  Five clusters: the weaker interval, printed as a
    twin to the i.i.d. one, never as a replacement for it."""
    d = np.asarray(d, float); folds = np.asarray(folds)
    fs = np.unique(folds)
    g = [d[folds == f] for f in fs]
    rng = np.random.default_rng(seed)
    bs = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, len(fs), len(fs))
        bs[b] = np.concatenate([g[i] for i in idx]).mean()
    return [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]


def pstat(a, b, folds, names=None):
    """Paired a - b with everything the brief requires attached."""
    p = I.paired(np.asarray(a, float), np.asarray(b, float), folds=folds, names=names)
    d = np.asarray(a, float) - np.asarray(b, float)
    p["ci95_fold"] = fold_boot(d, folds)
    p["conc"] = conc_verdict(d) if len(d) > 20 else None
    return p


def fmt(p):
    c = p["conc"]["verdict"] if p.get("conc") else "n/a"
    return (f"{p['mean_diff']:+.4f} [{p['ci95'][0]:+.4f}, {p['ci95'][1]:+.4f}] "
            f"fold[{p['ci95_fold'][0]:+.4f}, {p['ci95_fold'][1]:+.4f}] "
            f"med {p['median_diff']:+.4f} {p['n_better']}W/{p['n_worse']}L "
            f"n={p['n']} conc {c}")


# ------------------------------------------------------------------ Q1
def q1(rows, out, setnames):
    SETS = [BY_NAME[nm] for nm in setnames]
    folds = np.array([r["fold"] for r in rows])
    names = [r["pdb"] for r in rows]
    proj = np.array([r["proj"]["rmsd"] for r in rows])
    raw = np.array([r["raw_avg_rmsd"] for r in rows])
    n = len(rows)

    print("=" * 100)
    print("Q1.  RESTRAINED AMBER MINIMISATION AS A REPAIR OPERATOR -- the (k, depth) sweep")
    print(f"     n = {n} targets, TARGET is the unit, convergence gate = final energy finite "
          f"and <= {CONVERGE_MAX_KCAL:g} kcal/mol")
    print("=" * 100)

    # ---- the reference points, none of which is an AMBER arm
    print("\nREFERENCE POINTS (no AMBER):")
    print(f"  ZERO-INFORMATION 'do nothing'  raw all-atom coordinate average   {raw.mean():.4f} A"
          f"   median {np.median(raw):.4f}")
    print(f"  incumbent repair operator      ideal-geometry projection          {proj.mean():.4f} A"
          f"   median {np.median(proj):.4f}")
    hel = np.array([r["helix"]["rmsd"] for r in rows])
    print(f"  ZERO-INFORMATION constant ideal alpha-helix                       {hel.mean():.4f} A"
          f"   median {np.median(hel):.4f}")
    p = pstat(proj, raw, folds, names)
    print(f"  projection - do nothing : {fmt(p)}")
    out["ref"] = dict(raw=float(raw.mean()), proj=float(proj.mean()), helix=float(hel.mean()),
                      proj_minus_raw=p)

    hdr = (f"\n{'setting':>12s} {'k':>5s} {'depth':>6s} | {'AMBER':>7s} {'excl':>4s} | "
           f"{'vs PROJECTION (accuracy)':>46s} | {'vs DO-NOTHING':>22s} | "
           f"{'FRAME NULL (zero by construction)':>34s} | {'stable-subset effect':>26s}")
    print(hdr)
    out["settings"] = {}
    for nm, k, steps in SETS:
        A0 = np.array([r[f"f0_{nm}"]["rmsd"] for r in rows])
        A1 = np.array([r[f"f1_{nm}"]["rmsd"] for r in rows])
        E0 = np.array([r[f"f0_{nm}"]["energy"] for r in rows])
        E1 = np.array([r[f"f1_{nm}"]["energy"] for r in rows])
        c0 = np.array([r[f"f0_{nm}"]["converged"] for r in rows])
        c1 = np.array([r[f"f1_{nm}"]["converged"] for r in rows])
        keep = c0                                  # the gate, on the reported (frame-0) arm
        both = c0 & c1
        stable = both & (np.abs(E1 - E0) <= STABLE_DE)

        rec = dict(k=k, steps=steps, mean_rmsd=float(A0.mean()),
                   median_rmsd=float(np.median(A0)),
                   n_excluded=int((~keep).sum()),
                   excluded=[names[i] for i in np.where(~keep)[0]],
                   wall_mean=float(np.mean([r[f"f0_{nm}"]["wall"] for r in rows])),
                   moved_ca=float(np.mean([r[f"f0_{nm}"]["moved_ca"] for r in rows])))
        rec["vs_proj_ungated"] = pstat(A0, proj, folds, names)
        rec["vs_raw_ungated"] = pstat(A0, raw, folds, names)
        if keep.sum() >= 10:
            rec["vs_proj_gated"] = pstat(A0[keep], proj[keep], folds[keep],
                                         [names[i] for i in np.where(keep)[0]])
            rec["vs_raw_gated"] = pstat(A0[keep], raw[keep], folds[keep],
                                        [names[i] for i in np.where(keep)[0]])
        # the exact null
        d = A1 - A0
        nrec = {}
        for tag, m in (("ungated", np.ones(n, bool)), ("gated", both)):
            if m.sum() < 3:
                continue
            nrec[tag] = dict(n=int(m.sum()), mean=float(d[m].mean()),
                             sd=float(d[m].std(ddof=1)), max_abs=float(np.abs(d[m]).max()),
                             ci95=I.paired(A1[m], A0[m], folds=folds[m])["ci95"],
                             n_gt_1e6=int((np.abs(d[m]) > 1e-6).sum()))
        v = nrec.get("gated", nrec.get("ungated"))
        nrec["PASS"] = bool(abs(v["mean"]) <= FRAME_TOL_MEAN and v["max_abs"] <= FRAME_TOL_MAX)
        nrec["n_excluded_by_gate"] = int((~both).sum())
        rec["frame_null"] = nrec
        rec["n_stable"] = int(stable.sum())
        if stable.sum() >= 10:
            rec["vs_proj_stable"] = pstat(A0[stable], proj[stable], folds[stable],
                                          [names[i] for i in np.where(stable)[0]])
        out["settings"][nm] = rec

        vp = rec.get("vs_proj_gated", rec["vs_proj_ungated"])
        vr = rec.get("vs_raw_gated", rec["vs_raw_ungated"])
        vs = rec.get("vs_proj_stable")
        print(f"{nm:>12s} {k:5.0f} {('full' if steps == 0 else steps):>6} | "
              f"{A0.mean():7.4f} {rec['n_excluded']:4d} | "
              f"{vp['mean_diff']:+.4f} [{vp['ci95'][0]:+.4f},{vp['ci95'][1]:+.4f}] "
              f"{vp['n_better']:3d}W/{vp['n_worse']:3d}L | "
              f"{vr['mean_diff']:+.4f} {vr['n_better']:3d}W/{vr['n_worse']:3d}L | "
              f"mean {v['mean']:+.5f} max {v['max_abs']:.4f} n>0 {v['n_gt_1e6']:3d} "
              f"{'PASS' if nrec['PASS'] else 'FAIL'} | "
              + (f"{vs['mean_diff']:+.4f} [{vs['ci95'][0]:+.4f},{vs['ci95'][1]:+.4f}] "
                 f"{vs['n_better']:3d}W/{vs['n_worse']:3d}L n={rec['n_stable']}"
                 if vs else f"n={rec['n_stable']} too small"))

    # ---- the RANDOM-DIRECTION control at matched magnitude
    print("\nRANDOM-DIRECTION CONTROL at MATCHED MAGNITUDE "
          f"({N_RAND} draws, stable_rng; the input CA displaced by an isotropic random "
          "direction scaled to the SAME RMS displacement the operator applied)")
    print(f"{'setting':>12s} {'moved':>7s} | {'operator':>8s} {'random':>8s} "
          f"{'operator - random':>34s}")
    out["rand_control"] = {}
    for nm, k, steps in list(SETS) + [("proj", None, None)]:
        rr, oo = [], []
        for r in rows:
            A = np.array(r["input_ca"], float)
            nat = np.array(r["nat_ca"], float)
            nres = len(A)
            key = "proj" if nm == "proj" else f"f0_{nm}"
            mv = r[key]["moved_ca"]
            oo.append(r[key]["rmsd"])
            rng = stable_rng("s16repair", "randdir", nm, r["pdb"])
            vals = []
            for _ in range(N_RAND):
                g = rng.normal(size=(nres, 3))
                Q = AL.rigid_basis(A)
                g = (g.ravel() - Q @ (Q.T @ g.ravel())).reshape(nres, 3)
                g *= mv * np.sqrt(nres) / max(np.linalg.norm(g), 1e-12)
                vals.append(I.ca_rmsd(A + g, nat))
            rr.append(float(np.mean(vals)))
        rr = np.array(rr); oo = np.array(oo)
        p = pstat(oo, rr, folds, names)
        out["rand_control"][nm] = dict(operator=float(oo.mean()), random=float(rr.mean()),
                                       paired=p)
        mv = float(np.mean([r[("proj" if nm == "proj" else f"f0_{nm}")]["moved_ca"]
                            for r in rows]))
        print(f"{nm:>12s} {mv:7.4f} | {oo.mean():8.4f} {rr.mean():8.4f} "
              f"{p['mean_diff']:+.4f} [{p['ci95'][0]:+.4f},{p['ci95'][1]:+.4f}] "
              f"{p['n_better']:3d}W/{p['n_worse']:3d}L")

    # ---- VALIDITY, a SEPARATE axis
    print("\nVALIDITY -- a SEPARATE axis, never substituted for accuracy.  "
          "Rama = favoured fraction (higher better); clash = heavy atoms < 2.0 A; "
          "minhv = min heavy separation; geom = rms relative bond/angle deviation")
    print(f"{'arm':>14s} {'RMSD':>7s} {'rama':>6s} {'clash':>6s} {'minhv':>6s} {'geom':>7s} "
          f"{'wall s':>7s}")

    def vline(tag, rmsd, rama, geom, wall=None):
        print(f"{tag:>14s} {np.mean(rmsd):7.4f} {np.mean(rama):6.3f} "
              f"{np.mean([g['n_clash_2A'] for g in geom]):6.3f} "
              f"{np.mean([g['min_heavy'] for g in geom]):6.3f} "
              f"{np.mean([g['geom_rms_rel_dev'] for g in geom]):7.4f} "
              f"{('%7.2f' % np.mean(wall)) if wall is not None else '      -'}")

    out["validity"] = {}
    vline("do-nothing", raw, [r["raw_avg_rama_ok"] for r in rows],
          [r["raw_avg_geom"] for r in rows])
    vline("ZERO-INFO helix", hel, [r["helix"]["rama_ok"] for r in rows],
          [r["helix"]["geom"] for r in rows])
    vline("projection", proj, [r["proj"]["rama_ok"] for r in rows],
          [r["proj"]["geom"] for r in rows], [r["proj"]["wall"] for r in rows])
    for nm, k, steps in SETS:
        vline(nm, [r[f"f0_{nm}"]["rmsd"] for r in rows],
              [r[f"f0_{nm}"]["rama_ok"] for r in rows],
              [r[f"f0_{nm}"]["geom"] for r in rows],
              [r[f"f0_{nm}"]["wall"] for r in rows])
    for tag, get in (("raw_avg", lambda r: (r["raw_avg_rama_ok"], r["raw_avg_geom"])),
                     ("helix", lambda r: (r["helix"]["rama_ok"], r["helix"]["geom"])),
                     ("proj", lambda r: (r["proj"]["rama_ok"], r["proj"]["geom"]))):
        out["validity"][tag] = dict(
            rama=float(np.mean([get(r)[0] for r in rows])),
            clash=float(np.mean([get(r)[1]["n_clash_2A"] for r in rows])),
            min_heavy=float(np.mean([get(r)[1]["min_heavy"] for r in rows])),
            geom=float(np.mean([get(r)[1]["geom_rms_rel_dev"] for r in rows])))
    for nm, k, steps in SETS:
        out["validity"][nm] = dict(
            rama=float(np.mean([r[f"f0_{nm}"]["rama_ok"] for r in rows])),
            clash=float(np.mean([r[f"f0_{nm}"]["geom"]["n_clash_2A"] for r in rows])),
            min_heavy=float(np.mean([r[f"f0_{nm}"]["geom"]["min_heavy"] for r in rows])),
            geom=float(np.mean([r[f"f0_{nm}"]["geom"]["geom_rms_rel_dev"] for r in rows])),
            wall=float(np.mean([r[f"f0_{nm}"]["wall"] for r in rows])))
    # the paired validity comparisons that matter, against BOTH controls
    print("\nPAIRED VALIDITY (negative = the arm is better; Rama sign flipped so negative "
          "= better on every row)")
    for nm in [x for x in ("k30_full", "k0_full", "k60_full") if x in setnames]:
        for base_tag, base_r, base_g in (
                ("projection", [r["proj"]["rama_ok"] for r in rows],
                 [r["proj"]["geom"] for r in rows]),
                ("ZERO-INFO helix", [r["helix"]["rama_ok"] for r in rows],
                 [r["helix"]["geom"] for r in rows])):
            pr = pstat([-r[f"f0_{nm}"]["rama_ok"] for r in rows], [-x for x in base_r],
                       folds, names)
            pc = pstat([r[f"f0_{nm}"]["geom"]["n_clash_2A"] for r in rows],
                       [g["n_clash_2A"] for g in base_g], folds, names)
            print(f"  {nm} - {base_tag:>16s}  rama {pr['mean_diff']:+.4f} "
                  f"[{pr['ci95'][0]:+.4f},{pr['ci95'][1]:+.4f}] "
                  f"{pr['n_better']:3d}W/{pr['n_worse']:3d}L   "
                  f"clash {pc['mean_diff']:+.4f} "
                  f"[{pc['ci95'][0]:+.4f},{pc['ci95'][1]:+.4f}] "
                  f"{pc['n_better']:3d}W/{pc['n_worse']:3d}L")
            out["validity"].setdefault("paired", {})[f"{nm}_vs_{base_tag}"] = dict(
                rama=pr, clash=pc)
    return out


# ------------------------------------------------------------------ Q2
def _shares(x, Urot, h, rank, Q):
    """Fractions of ||x||^2 in the loud half / quiet half / non-torsional / rigid parts."""
    x = np.asarray(x, float)
    nx2 = float(x @ x)
    if nx2 <= 0:
        return dict(loud=np.nan, quiet=np.nan, nontors=np.nan, rigid=np.nan)
    c = Urot.T @ x
    loud = float((c[:h] ** 2).sum() / nx2)
    quiet = float((c[h:rank] ** 2).sum() / nx2)
    tot = float((c[:rank] ** 2).sum() / nx2)
    rig = float(((Q.T @ x) ** 2).sum() / nx2)
    return dict(loud=loud, quiet=quiet, nontors=float(1.0 - tot), rigid=rig)


def q2(rows, out, setnames):
    print("\n" + "=" * 100)
    print("Q2.  WHAT THE REPAIR OPERATORS DO TO THE GEOMETRY")
    print("     Every cosine and every residual below reads the NATIVE and is an ORACLE")
    print("     DIAGNOSTIC.  None of it enters any parameter, threshold or stopping rule.")
    print("=" * 100)
    ops = [x for x in ("k30_full", "k0_full", "k60_full", "k30_d25", "k30_d100")
           if x in setnames] + ["proj", "helix"] + [f"rand{i}" for i in range(N_RAND)]
    acc = {o: {kk: [] for kk in
               ("cos", "mag", "loud", "quiet", "nontors", "rigid", "gain_opt",
                "cos_tors", "loud_tors", "quiet_tors", "null_tors", "mag_tors")}
           for o in ops}
    res = {kk: [] for kk in ("loud", "quiet", "nontors", "rigid", "rmsd")}
    diffcos, diffmag, pairc = [], [], []
    # the ANALYTIC isotropic-direction null for subspace occupancy.  A random direction in
    # the (3n-6)-dimensional rigid-orthogonal CA space puts share dim/(3n-6) in each block;
    # in torsion space the denominators are 2n.  This is the RANDOM-DIRECTION control for
    # "which half of the spectrum does the operator move in", and it is exact, not sampled.
    nullsh = {kk: [] for kk in ("loud", "quiet", "nontors", "loud_t", "quiet_t", "null_t")}
    csteer_proj = []
    tg = {t["pdb"]: t for t in I.targets()}

    for r in rows:
        p = r["pdb"]
        W, PHI, PSI, u = top75_windows(p)
        avg, C_ca, _ = averaged_backbone_from(W, PHI, PSI)
        A = np.asarray(r["input_ca"], float)
        nat = np.asarray(r["nat_ca"], float)
        n = len(A)
        phi_in, psi_in = torsions_of(avg)
        nt = PDB.by_pdb(p)                                     # ORACLE torsions
        nphi, npsi = np.asarray(nt.phi, float), np.asarray(nt.psi, float)
        e_tors = np.concatenate([AL.wrap(nphi - phi_in), AL.wrap(npsi - psi_in)])

        J, _, CAJ = AL.sup_jacobian(phi_in, psi_in)
        U, s, Vt = np.linalg.svd(J, full_matrices=False)
        V = Vt.T
        rank = int((s > 1e-8 * s[0]).sum())
        h = rank // 2
        R, Am, Bm = _kabsch(CAJ[None], A[None])
        R = R[0]
        Urot = (U.T.reshape(-1, n, 3) @ R.T).reshape(-1, 3 * n).T
        Q = AL.rigid_basis(A)

        natA = _sup(nat, A)
        rres = (natA - A).ravel()
        nr = float(np.linalg.norm(rres))
        sh = _shares(rres, Urot, h, rank, Q)
        for kk in ("loud", "quiet", "nontors", "rigid"):
            res[kk].append(sh[kk])
        res["rmsd"].append(nr / np.sqrt(n))

        d_ca = 3 * n - 6
        nullsh["loud"].append(h / d_ca)
        nullsh["quiet"].append((rank - h) / d_ca)
        nullsh["nontors"].append((d_ca - rank) / d_ca)
        nullsh["loud_t"].append(h / (2 * n))
        nullsh["quiet_t"].append((rank - h) / (2 * n))
        nullsh["null_t"].append((2 * n - rank) / (2 * n))

        # CROSS-CHECK against `s16/csteer.py`, which measured the same projection direction
        # from a DIFFERENT start object: its `avg` is `I.coordinate_average(W)` (CA only),
        # not the CA block of the all-atom average this operator is handed.  Both are
        # computed here so the two workstreams' numbers can be reconciled rather than
        # silently disagreeing.
        Cs = np.asarray(C_ca, float)
        rc = (_sup(nat, Cs) - Cs).ravel()
        vc = (_sup(np.asarray(r["proj"]["ca"], float), Cs) - Cs).ravel()
        csteer_proj.append(float(vc @ rc / max(np.linalg.norm(vc) * np.linalg.norm(rc), 1e-12)))

        rng = stable_rng("s16repair", "q2rand", p)
        vstore = {}
        for o in ops:
            if o.startswith("rand"):
                g = rng.normal(size=3 * n)
                g = g - Q @ (Q.T @ g)
                ref = np.asarray(r["f0_k30_full"]["ca"], float)
                mv = float(np.linalg.norm((_sup(ref, A) - A).ravel()))
                v = g * mv / max(np.linalg.norm(g), 1e-12)
                dth = None
            else:
                key = o if o in ("proj", "helix") else f"f0_{o}"
                B = np.asarray(r[key]["ca"], float)
                v = (_sup(B, A) - A).ravel()
                if o == "helix":
                    dth = np.concatenate([AL.wrap(np.full(n, np.radians(-57.0)) - phi_in),
                                          AL.wrap(np.full(n, np.radians(-47.0)) - psi_in)])
                elif "phi" in r[key]:
                    dth = np.concatenate([AL.wrap(np.asarray(r[key]["phi"], float) - phi_in),
                                          AL.wrap(np.asarray(r[key]["psi"], float) - psi_in)])
                else:
                    dth = None
            vstore[o] = v
            nv = float(np.linalg.norm(v))
            c = float(v @ rres / max(nv * nr, 1e-12))
            acc[o]["cos"].append(c)
            acc[o]["mag"].append(nv / max(nr, 1e-12))
            acc[o]["gain_opt"].append(nr * (1.0 - np.sqrt(max(0.0, 1.0 - c * c))) / np.sqrt(n))
            shv = _shares(v, Urot, h, rank, Q)
            for kk in ("loud", "quiet", "nontors", "rigid"):
                acc[o][kk].append(shv[kk])
            if dth is not None:
                nd = float(np.linalg.norm(dth)); ne = float(np.linalg.norm(e_tors))
                acc[o]["cos_tors"].append(float(dth @ e_tors / max(nd * ne, 1e-12)))
                acc[o]["mag_tors"].append(nd / max(ne, 1e-12))
                cc = V.T @ dth
                q2s = float(cc @ cc)
                acc[o]["loud_tors"].append(float((cc[:h] ** 2).sum() / q2s))
                acc[o]["quiet_tors"].append(float((cc[h:rank] ** 2).sum() / q2s))
                acc[o]["null_tors"].append(float((cc[rank:] ** 2).sum() / q2s))

        # the DIFFERENCE between the two repair operators -- this, not either displacement,
        # is the direction the programme's quoted -0.022 A effect actually lives along.
        if "k30_full" in vstore and "proj" in vstore:
            dv = vstore["k30_full"] - vstore["proj"]
            nd = float(np.linalg.norm(dv))
            diffcos.append(float(dv @ rres / max(nd * nr, 1e-12)))
            diffmag.append(nd / np.sqrt(n))
            pairc.append(float(vstore["k30_full"] @ vstore["proj"] /
                               max(np.linalg.norm(vstore["k30_full"]) *
                                   np.linalg.norm(vstore["proj"]), 1e-12)))

    folds = np.array([r["fold"] for r in rows])
    names = [r["pdb"] for r in rows]
    print("\nA.  THE TRUE RESIDUAL (ORACLE) at the input structure, and where it lives")
    print(f"    ||r||/sqrt(n) = {np.mean(res['rmsd']):.4f} A   "
          f"loud half {np.mean(res['loud']):.3f}  quiet half {np.mean(res['quiet']):.3f}  "
          f"non-torsional {np.mean(res['nontors']):.3f}  rigid {np.mean(res['rigid']):.4f}")
    out["q2_residual"] = {k: float(np.mean(v)) for k, v in res.items()}

    print("\nB.  THE DISPLACEMENT EACH OPERATOR APPLIES, in CA space "
          "(superposed onto the input; rigid part removed by the superposition)")
    print(f"{'operator':>10s} {'cos(v,r)':>22s} {'|v|/|r|':>8s} {'loud':>6s} {'quiet':>6s} "
          f"{'nontor':>7s} | {'gain if optimally stepped':>26s}")
    out["q2_ops"] = {}
    for o in ops:
        c = np.array(acc[o]["cos"])
        ci = I.paired(c, np.zeros_like(c), folds=folds)["ci95"]
        g = np.array(acc[o]["gain_opt"])
        out["q2_ops"][o] = {k: float(np.mean(v)) for k, v in acc[o].items() if len(v)}
        out["q2_ops"][o]["cos_ci95"] = [float(x) for x in ci]
        out["q2_ops"][o]["cos_median"] = float(np.median(c))
        out["q2_ops"][o]["cos_frac_pos"] = float((c > 0).mean())
        print(f"{o:>10s} {c.mean():+.4f} [{ci[0]:+.4f},{ci[1]:+.4f}] med {np.median(c):+.3f} "
              f"{np.mean(acc[o]['mag']):8.4f} {np.mean(acc[o]['loud']):6.3f} "
              f"{np.mean(acc[o]['quiet']):6.3f} {np.mean(acc[o]['nontors']):7.3f} | "
              f"{np.mean(g):.4f} A  (frac cos>0 {(c > 0).mean():.3f})")
    print(f"{'ANALYTIC':>10s} {'(isotropic-direction null)':>40s} {'-':>8s} "
          f"{np.mean(nullsh['loud']):6.3f} {np.mean(nullsh['quiet']):6.3f} "
          f"{np.mean(nullsh['nontors']):7.3f} |  <- the RANDOM-DIRECTION control for "
          f"subspace occupancy")
    out["q2_subspace_null_ca"] = {k: float(np.mean(nullsh[k]))
                                  for k in ("loud", "quiet", "nontors")}

    print("\nC.  THE SAME DISPLACEMENT IN TORSION SPACE, decomposed in the RIGHT singular "
          "vectors of the superposed-CA Jacobian (loud first)")
    print(f"{'operator':>10s} {'cos(dtheta, e_true)':>28s} {'|dth|/|e|':>10s} "
          f"{'loud':>6s} {'quiet':>6s} {'exact-null':>10s}")
    for o in ops:
        if not acc[o]["cos_tors"]:
            continue
        c = np.array(acc[o]["cos_tors"])
        ci = I.paired(c, np.zeros_like(c), folds=folds)["ci95"]
        print(f"{o:>10s} {c.mean():+.4f} [{ci[0]:+.4f},{ci[1]:+.4f}] med {np.median(c):+.3f} "
              f"{np.mean(acc[o]['mag_tors']):10.4f} {np.mean(acc[o]['loud_tors']):6.3f} "
              f"{np.mean(acc[o]['quiet_tors']):6.3f} {np.mean(acc[o]['null_tors']):10.3f}")
    print(f"{'ANALYTIC':>10s} {'(isotropic-direction null)':>28s} {'-':>10s} "
          f"{np.mean(nullsh['loud_t']):6.3f} {np.mean(nullsh['quiet_t']):6.3f} "
          f"{np.mean(nullsh['null_t']):10.3f}")
    out["q2_subspace_null_tors"] = {k: float(np.mean(nullsh[k]))
                                    for k in ("loud_t", "quiet_t", "null_t")}
    cp = np.array(csteer_proj)
    cip = I.paired(cp, np.zeros_like(cp), folds=folds)["ci95"]
    print("")
    print("    CROSS-CHECK vs s16/csteer.py, which measured cos(avg -> proj, r) = -0.0677")
    print("    from the CA-only coordinate average.  Recomputed here from THAT start: "
          f"{cp.mean():+.4f} [{cip[0]:+.4f}, {cip[1]:+.4f}], n = {len(cp)}.")
    print("    From the ALL-ATOM average's CA block -- the object AMBER is actually handed,")
    print("    which differs from it by ~0.16 A -- the same direction reads "
          f"{np.mean(acc['proj']['cos']):+.4f}.")
    out["q2_csteer_crosscheck"] = dict(from_C_ca=float(cp.mean()),
                                       ci95=[float(x) for x in cip],
                                       csteer_reported=-0.0677,
                                       from_allatom_avg=float(np.mean(acc['proj']['cos'])))

    if diffcos:
        dc = np.array(diffcos)
        ci = I.paired(dc, np.zeros_like(dc), folds=folds)["ci95"]
        print("\nF.  THE DIRECTION THE QUOTED EFFECT LIVES ALONG.  The -0.022 A figure is a")
        print("    DIFFERENCE of two displacements, v_AMBER - v_projection.  If that difference")
        print("    has no cosine with the true residual, the effect is not physics finding a")
        print("    better structure -- it is one operator paying less for valid geometry.")
        print(f"    cos(v_AMBER - v_proj, r) = {dc.mean():+.4f} [{ci[0]:+.4f}, {ci[1]:+.4f}]  "
              f"median {np.median(dc):+.4f}  frac>0 {(dc > 0).mean():.3f}")
        print(f"    ||v_AMBER - v_proj||/sqrt(n) = {np.mean(diffmag):.4f} A   "
              f"cos(v_AMBER, v_proj) = {np.mean(pairc):+.4f}")
        out["q2_diff"] = dict(cos=float(dc.mean()), cos_ci95=[float(x) for x in ci],
                              cos_median=float(np.median(dc)),
                              frac_pos=float((dc > 0).mean()),
                              mag=float(np.mean(diffmag)), cos_ops=float(np.mean(pairc)))
        out["q2_ops"][o]["cos_tors_ci95"] = [float(x) for x in ci]

    print("\nE.  IS THE OPERATOR'S RMSD CHANGE EXPLAINED BY ITS DISPLACEMENT ALONE?")
    print("    In coordinate space  n*RMSD_after^2 = ||r||^2 - 2 v.r + ||v||^2  exactly, so a")
    print("    displacement ORTHOGONAL to the residual can only cost, and the cost is fixed by")
    print("    its magnitude: dRMSD ~ ||v||^2 / (2 n RMSD).  Predicted vs realised, n = %d."
          % len(rows))
    print(f"{'operator':>10s} {'RMSD in':>8s} {'predicted':>10s} {'realised':>9s} "
          f"{'|pred-real|':>11s} {'cost of |v| alone':>18s} {'corr(|v|, dRMSD)':>17s}")
    nres = np.array([len(r["input_ca"]) for r in rows])
    r0 = np.array(res["rmsd"])
    for o in ops:
        if o.startswith("rand"):
            continue
        key = o if o in ("proj", "helix") else f"f0_{o}"
        real = np.array([r[key]["rmsd"] for r in rows])
        mag = np.array(acc[o]["mag"]) * r0 * np.sqrt(nres)          # ||v||
        cs = np.array(acc[o]["cos"])
        pred = np.sqrt(np.maximum(0.0, ((r0 * np.sqrt(nres)) ** 2
                                        - 2 * cs * mag * r0 * np.sqrt(nres)
                                        + mag ** 2)) / nres)
        cost = (mag ** 2) / (2 * nres * r0)                          # the orthogonal-move cost
        cc = np.corrcoef(mag / np.sqrt(nres), real - r0)[0, 1]
        out["q2_ops"][o]["pred_rmsd"] = float(pred.mean())
        out["q2_ops"][o]["real_rmsd"] = float(real.mean())
        out["q2_ops"][o]["orth_cost"] = float(cost.mean())
        out["q2_ops"][o]["corr_mag_drmsd"] = float(cc)
        print(f"{o:>10s} {r0.mean():8.4f} {pred.mean():10.4f} {real.mean():9.4f} "
              f"{np.abs(pred - real).mean():11.4f} {cost.mean():18.4f} {cc:17.3f}")

    print("\nD.  PAIRED, against the two mandatory controls (positive = the operator's "
          "direction is better aligned with the truth)")
    base = {"ZERO-INFO helix": "helix", "RANDOM matched-magnitude": "rand0"}
    for o in [x for x in ("k30_full", "k0_full", "proj") if x in ops]:
        for bt, bo in base.items():
            p = pstat(acc[o]["cos"], acc[bo]["cos"], folds, names)
            print(f"  cos(v,r):  {o:>10s} - {bt:<26s} {fmt(p)}")
            out.setdefault("q2_paired", {})[f"{o}_vs_{bo}"] = p
    return out


def main(passes=("A", "B", "C")):
    setnames, rows = load(passes)
    out = dict(what="REPAIR: AMBER as an operator", n=len(rows), passes=list(passes),
               settings=setnames)
    q1(rows, out, setnames)
    q2(rows, out, setnames)
    with open(os.path.join(RESULTS, "repair_report.json"), "w") as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("\nwrote s16/results/repair_report.json")


if __name__ == "__main__":
    ps = sys.argv[1].split(",") if len(sys.argv) > 1 else ("A", "B", "C")
    main(tuple(ps))
