"""SPRINT 14, ENER-10 -- can PHYSICS replace the ideal-geometry PROJECTION?

The coordinator measured (`s14/avgspace.py`) that on the 126-target instrument, with the
shipped top-75 windows held constant in every arm:

    A  coordinate average -> project     3.205
    B  coordinate average, raw           3.048
    projection costs +0.157 A [+0.124, +0.191], losing on 108 of 126 targets.

The proposal: the projection buys validity by snapping to IDEAL geometry (fixed bonds,
fixed angles, only torsions free).  A force field needs no such thing.  Relax arm B with
restrained AMBER instead and see whether validity can be recovered for less than 0.157 A.

## The confound in the test as posed, and the fix

`I.coordinate_average` averages `u["W"]`, which is a **CA trace only** -- (n, 3), no
backbone.  `core.amber.refine_coords` requires N, CA, C, O and CB.  The only route in this
repository from a bare CA trace to an all-atom structure is `core.amber.refine_ca`, and
that function **fits the discrete torsion states to the CA trace first** (`floor.descend`).
That fit IS a projection onto the ideal-geometry manifold -- the exact step under test.
Running the experiment that way would measure projection-then-AMBER and report it as
AMBER-instead-of-projection.  **The test as literally posed is not runnable.**

The fix keeps the experiment intact and makes it stronger.  Every window carries its own
`PHI`/`PSI`, so each window has a full ideal-geometry backbone.  Superpose each window's
full backbone onto that window's own CA trace IN THE MEDOID FRAME -- the same frame
`coordinate_average` works in -- then average every atom.  The result is a genuine averaged
full backbone whose bond lengths and angles are whatever averaging produced (non-ideal,
exactly the object the proposal wants) and whose CA block reproduces arm B.  The agreement
is ASSERTED AT RUN TIME and reported, never assumed: `ca_identity_max_dev` is the max
coordinate deviation and the two arms' RMSDs are printed side by side.  No projection
anywhere.

## THE PREDICTION, stated before the measurement

Averaging superposed structures CONTRACTS them: the mean of a set of unit vectors is
shorter than a unit vector, so the averaged backbone's bonds must come out short, and the
more the windows disagree the shorter they get.  If that is the mechanism, then the
+0.157 A is the price of RE-EXPANDING a shrunken structure, any operation that restores
valid bond lengths must pay it, and AMBER will land near 3.2 rather than near 3.05.
Outcome 2.  The contraction is measured first (`bond_contraction`) and its correlation with
the per-target projection cost is reported, so the mechanism is falsifiable independently
of the RMSD result.

## ARMS, all on the same 126 targets and the same windows

    raw_avg           arm B, the averaged full backbone, unrefined
    projected         arm A, the incumbent ideal-geometry projection
    legacy_refine     coordinate descent on the Legacy total from the projection's torsions
    amber_k10         restrained AMBER relaxation, k = 10 (the production restraint)
    amber_k2          restrained AMBER relaxation, k = 2  (loose)
    amber_k0          FREE AMBER relaxation, no restraint at all

The restraint strength is the knob that trades RMSD against validity, so reporting a single
k would be picking a point on a frontier and calling it a result.  All three are reported
with the geometric audit beside them, because an RMSD gain bought with invalid geometry is
not a gain.

    python -m s14.ener_avgrefine [--pilot N] [--k 10,2,0]
"""
from __future__ import annotations

import sys
import time

import numpy as np

from s12 import instrument as I
from s14 import ener_lib as E
from s14.avgspace import top75_windows

IDEAL = dict(N_CA=1.458, CA_C=1.525, C_N=1.329, C_O=1.231,
             ang_N_CA_C=111.0, ang_CA_C_N=116.2, ang_C_N_CA=121.7, ca_ca=3.804)
ATOMS = ("N", "CA", "C", "O", "CB")


# ------------------------------------------------------------------ arm B, full backbone
def _kabsch(A, B):
    """Rigid transform taking each A[b] onto B[b]. Returns (R, tA, tB)."""
    Am = A.mean(1, keepdims=True); Bm = B.mean(1, keepdims=True)
    Ac = A - Am; Bc = B - Bm
    H = np.einsum("bni,bnj->bij", Ac, Bc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("bij,bjk->bik", Vt.transpose(0, 2, 1),
                                        U.transpose(0, 2, 1))))
    D = np.tile(np.eye(3), (len(A), 1, 1)); D[:, 2, 2] = d
    R = np.einsum("bij,bjk,bkl->bil", Vt.transpose(0, 2, 1), D, U.transpose(0, 2, 1))
    return R, Am, Bm


def averaged_backbone(pdb):
    """Arm B extended to every backbone atom.

    CORRECTION to the first version of this function, which was wrong and whose own
    identity assertion caught it (16.1 A deviation).  `geo.build_backbone_batch(PHI, PSI)`
    returns each window's IDEAL-GEOMETRY REBUILD in the builder's own reference frame; it
    is not the window's stored CA coordinates `W` and it is not in `W`'s frame.  Applying
    the W-derived superposition to it superposes nothing and averages structures in
    unrelated frames, which is what produced an apparent 86% bond contraction.  That
    number was an artefact of the bug and is retracted.

    The correct construction, in two stages:

      1. `Wm = superpose_batch(W, W[medoid])` -- the real CA traces in the medoid frame.
         `Wm.mean(0)` is exactly arm B.
      2. Superpose each window's BUILT full backbone onto its OWN `Wm[i]` using the CA
         blocks, then average every atom.

    The CA block of the result differs from arm B only by the per-window ideal-geometry
    rebuild residual (the project's 0.347 A floor), most of which cancels in the mean.
    The actual deviation is returned as `ca_identity_max_dev` and reported, not assumed.
    """
    from core import geometry as geo
    W, PHI, PSI, u = top75_windows(pdb)
    C_ca, b = I.coordinate_average(W)               # the incumbent arm-B object
    Wm = I.superpose_batch(W, W[b])                 # real CAs in the medoid frame
    full = geo.build_backbone_batch(PHI, PSI)       # per-window ideal backbones
    bca = np.asarray(full["CA"], float)
    R, Am, Bm = _kabsch(bca, Wm)                    # built frame -> medoid frame
    avg = {}
    for a in ATOMS:
        if a not in full:
            continue
        X = np.asarray(full[a], float)
        avg[a] = (np.einsum("bij,bnj->bni", R, X - Am) + Bm).mean(0)
    dev = float(np.abs(avg["CA"] - C_ca).max())
    # NOT the reconstruction error -- this is the spread of the individual windows about
    # their own mean, reported because it is what the averaging has to cancel.
    spread = float(np.mean(I.kabsch_rmsd_batch(
        np.einsum("bij,bnj->bni", R, bca - Am) + Bm, C_ca)))
    return avg, C_ca, u, dev, len(W), spread


# ------------------------------------------------------------------ geometry
def geom_of(c, seq):
    from s14.ener_geom import _ang, _dih
    N, CA, Cc = (np.asarray(c[k], float) for k in ("N", "CA", "C"))
    O = np.asarray(c["O"], float) if "O" in c else None
    CB = np.asarray(c["CB"], float) if "CB" in c else None
    g = dict(N_CA=np.linalg.norm(CA - N, axis=1),
             CA_C=np.linalg.norm(Cc - CA, axis=1),
             C_N=np.linalg.norm(N[1:] - Cc[:-1], axis=1),
             ca_ca=np.linalg.norm(CA[1:] - CA[:-1], axis=1),
             ang_N_CA_C=_ang(N, CA, Cc),
             ang_CA_C_N=_ang(CA[:-1], Cc[:-1], N[1:]),
             ang_C_N_CA=_ang(Cc[:-1], N[1:], CA[1:]))
    if O is not None:
        g["C_O"] = np.linalg.norm(O - Cc, axis=1)
    out = {k: float(np.mean(v)) for k, v in g.items()}
    out.update({k + "_min": float(np.min(v)) for k, v in g.items()})
    out["omega"] = float(np.mean(np.abs(np.degrees(
        _dih(CA[:-1], Cc[:-1], N[1:], CA[1:])))))
    if CB is not None:
        ok = np.array([i for i, a in enumerate(seq) if a != "G"], int)
        if len(ok):
            imp = np.degrees(_dih(N[ok], CA[ok], Cc[ok], CB[ok]))
            out["chirality_mean"] = float(imp.mean())
            out["chirality_L_frac"] = float((imp < 0).mean())
    # non-bonded clash: heavy atoms at least two residues apart
    arrs, res = [], []
    for a in ATOMS:
        if a in c:
            v = np.asarray(c[a], float); arrs.append(v); res.append(np.arange(len(v)))
    Hh = np.vstack(arrs); Ri = np.concatenate(res)
    Dm = np.linalg.norm(Hh[:, None, :] - Hh[None, :, :], axis=-1)
    Dm = np.where(np.abs(Ri[:, None] - Ri[None, :]) >= 2, Dm, 9e9)
    out["min_heavy"] = float(Dm.min())
    out["n_clash_2A"] = int((Dm < 2.0).sum() // 2)
    # a single scalar: RMS relative deviation of every bond and angle from ideal
    dev = []
    for k, ideal in IDEAL.items():
        if k in out:
            dev.append((out[k] - ideal) / ideal)
    out["geom_rms_rel_dev"] = float(np.sqrt(np.mean(np.square(dev))))
    out["bond_contraction"] = float(np.mean([
        (out[k] - IDEAL[k]) / IDEAL[k] for k in ("N_CA", "CA_C", "C_N", "ca_ca")
        if k in out]))
    return out


def torsions_of(c):
    """phi/psi of an arbitrary backbone, so an averaged structure can be re-scored."""
    from s14.ener_geom import _dih
    N, CA, Cc = (np.asarray(c[k], float) for k in ("N", "CA", "C"))
    n = len(CA)
    phi = np.zeros(n); psi = np.zeros(n)
    phi[1:] = _dih(Cc[:-1], N[1:], CA[1:], Cc[1:])
    psi[:-1] = _dih(N[:-1], CA[:-1], Cc[:-1], N[1:])
    return phi, psi


def rama_ok(phi, psi):
    p = np.degrees(phi)[1:-1]; s = np.degrees(psi)[1:-1]
    a = (p > -160) & (p < -20) & (s > -120) & (s < 50)
    b = (p > -180) & (p < -40) & (s > 90) & (s < 180)
    l = (p > 20) & (p < 100) & (s > -20) & (s < 90)
    return float((a | b | l).mean()) if len(p) else float("nan")


# ------------------------------------------------------------------ the arms
def run_target(t, ks, seed=0):
    from core import amber as am
    import torsion_lib2 as tl2
    pdb, seq, fold = t["pdb"], t["seq"], int(t["fold"])
    avg, C_ca, u, dev, nwin, spread = averaged_backbone(pdb)
    nat = u["nat_ca"]
    rec = dict(pdb=pdb, fold=fold, n=len(seq), n_windows=int(nwin),
               ca_identity_max_dev=dev, window_spread=spread,
               armB_rmsd=float(I.ca_rmsd(C_ca, nat)),
               raw_avg_rmsd=float(I.ca_rmsd(avg["CA"], nat)),
               raw_avg_geom=geom_of(avg, seq))
    ph, ps = torsions_of(avg)
    rec["raw_avg_rama_ok"] = rama_ok(ph, ps)
    # arm A: the incumbent projection
    pr = I.project(C_ca, seq, fold)
    rec["projected_rmsd"] = float(I.ca_rmsd(pr["fit_ca"], nat))
    rec["projected_lam_rmsd"] = float(I.ca_rmsd(pr["ca"], nat))
    rec["projected_rama_ok"] = rama_ok(np.asarray(pr["phi"]), np.asarray(pr["psi"]))
    # a torsion representation for the AMBER builder (sidechain placement only; the
    # BACKBONE handed to AMBER is the averaged one, never this)
    tab = tl2.library_for(seq, 4, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    for k in ks:
        E.wait_for_memory(1.5, f"avgrefine-{pdb}")
        t0 = time.time()
        try:
            r = am.refine_coords(seq, rep, avg, k_restraint=float(k), steps=0,
                                 tolerance=1.0, threads=1, memo=False)
        except NotImplementedError as exc:      # a residue the sidechain builder lacks
            rec[f"amber_k{k}"] = dict(error=str(exc)[:120])
            continue
        ca = np.asarray(r["ca"], float)
        bb = {a: np.asarray(r["backbone"][a], float) for a in ATOMS
              if a in r["backbone"]}
        phi2, psi2 = torsions_of(bb)
        rec[f"amber_k{k}"] = dict(
            rmsd=float(I.ca_rmsd(ca, nat)), energy=float(r["energy"]),
            energy_initial=float(r["energy_initial"]),
            restraint_rmsd=float(r["restraint_rmsd"]),
            moved_ca=float(I.ca_rmsd(ca, avg["CA"])),
            rama_ok=rama_ok(phi2, psi2), geom=geom_of(bb, seq),
            wall=round(time.time() - t0, 2))
    return rec


# ------------------------------------------------------------------ main
def main(pilot=None, ks=(10.0, 2.0, 0.0)):
    tg = I.targets()
    if pilot:
        tg = tg[:int(pilot)]
    rows = []
    for i, t in enumerate(tg):
        rows.append(run_target(t, ks))
        r = rows[-1]
        msg = (f"[{i+1}/{len(tg)}] {r['pdb']} raw {r['raw_avg_rmsd']:.3f} "
               f"proj {r['projected_rmsd']:.3f}")
        for k in ks:
            a = r.get(f"amber_k{k}")
            if a and "rmsd" in a:
                msg += f" k{k:g} {a['rmsd']:.3f}(g{a['geom']['geom_rms_rel_dev']:.3f})"
        print(msg + f"  ident {r['ca_identity_max_dev']:.1e}", flush=True)
        if (i + 1) % 10 == 0 or i + 1 == len(tg):
            E.write("ener_avgrefine", dict(
                what="can restrained AMBER replace the ideal-geometry projection?",
                ks=list(ks), per_target=rows), n_expected=len(tg))
    report(rows, ks)
    return rows


def report(rows, ks):
    ok = [r for r in rows if all(f"amber_k{k}" in r and "rmsd" in r[f"amber_k{k}"]
                                 for k in ks)]
    folds = [r["fold"] for r in ok]
    f18 = set(I.FAIL18)
    print(f"\n=== ARM B + AMBER vs ARM A ({len(ok)} of {len(rows)} targets complete) ===")
    print(f"reconstructed arm B vs the true arm B: max coordinate deviation "
          f"{max(r['ca_identity_max_dev'] for r in rows):.3f} A, "
          f"mean window spread about the average {np.mean([r['window_spread'] for r in rows]):.3f} A")
    print(f"true arm B (CA average)  {np.mean([r['armB_rmsd'] for r in ok]):.3f} A   "
          f"reconstructed arm B (all-atom average) "
          f"{np.mean([r['raw_avg_rmsd'] for r in ok]):.3f} A")
    print(f"\n{'arm':16s} {'RMSD':>7s} {'FAIL18':>7s} {'other':>7s} {'vs A':>8s} "
          f"{'CI':>20s} {'W/L':>8s} {'geom dev':>9s} {'rama':>6s} {'clash':>6s}")
    A = np.array([r["projected_rmsd"] for r in ok])

    def line(name, v, geom, rama, clash):
        p = I.paired(np.asarray(v), A, folds=folds)
        fa = [x for x, r in zip(v, ok) if r["pdb"] in f18]
        ot = [x for x, r in zip(v, ok) if r["pdb"] not in f18]
        print(f"{name:16s} {np.mean(v):7.3f} {np.mean(fa):7.3f} {np.mean(ot):7.3f} "
              f"{p['mean_diff']:+8.3f} [{p['ci95'][0]:+.3f},{p['ci95'][1]:+.3f}] "
              f"{p['n_better']:3d}/{p['n_worse']:<3d} "
              f"{np.mean(geom):9.4f} {np.mean(rama):6.3f} {np.mean(clash):6.3f}")

    line("A projected", A, [0.0] * len(ok),
         [r["projected_rama_ok"] for r in ok], [0.0] * len(ok))
    line("B raw average", [r["raw_avg_rmsd"] for r in ok],
         [r["raw_avg_geom"]["geom_rms_rel_dev"] for r in ok],
         [r["raw_avg_rama_ok"] for r in ok],
         [float(r["raw_avg_geom"]["n_clash_2A"] > 0) for r in ok])
    for k in ks:
        line(f"B + AMBER k={k:g}", [r[f"amber_k{k}"]["rmsd"] for r in ok],
             [r[f"amber_k{k}"]["geom"]["geom_rms_rel_dev"] for r in ok],
             [r[f"amber_k{k}"]["rama_ok"] for r in ok],
             [float(r[f"amber_k{k}"]["geom"]["n_clash_2A"] > 0) for r in ok])
    print("\n--- THE PREDICTION: is the averaged backbone CONTRACTED? ---")
    for key in ("N_CA", "CA_C", "C_N", "ca_ca", "C_O"):
        v = [r["raw_avg_geom"][key] for r in ok if key in r["raw_avg_geom"]]
        if v:
            print(f"  {key:8s} averaged {np.mean(v):7.4f}  ideal {IDEAL[key]:7.4f}  "
                  f"{100*(np.mean(v)-IDEAL[key])/IDEAL[key]:+6.2f}%")
    bc = [r["raw_avg_geom"]["bond_contraction"] for r in ok]
    print(f"  mean bond contraction {100*np.mean(bc):+.2f}%   "
          f"per-target range [{100*min(bc):+.2f}%, {100*max(bc):+.2f}%]")
    print(f"  correlation( contraction , A-minus-B projection cost ) = "
          f"{E.spearman(bc, [r['projected_rmsd'] - r['raw_avg_rmsd'] for r in ok]):+.3f}")
    print("\n--- how far did AMBER move the structure? ---")
    for k in ks:
        mv = [r[f"amber_k{k}"]["moved_ca"] for r in ok]
        print(f"  k={k:<5g} CA moved {np.mean(mv):.3f} A   "
              f"wall {np.mean([r[f'amber_k{k}']['wall'] for r in ok]):.1f} s")
    E.write("ener_avgrefine", dict(
        what="can restrained AMBER replace the ideal-geometry projection?",
        ks=list(ks), per_target=rows), n_expected=len(rows))


def concentration_null(d, n_sim=4000, seed=0):
    """What the drop-top statistics look like under a UNIFORM effect of the same size.

    The drop-top curve is the project's mandated concentration check, but it is only
    interpretable against this null.  For an effect whose mean is small compared with the
    per-target spread, discarding the 10 most favourable targets removes a large share of
    the total even when every target carries an identical effect -- so a raw threshold
    reports a uniform effect as concentrated.  Simulating N(mean, sd) draws of the same
    length gives the reference distribution.
    """
    d = np.asarray(d, float)
    n = len(d)
    rng = np.random.default_rng(seed)
    sh, a10, a20 = [], [], []
    for _ in range(n_sim):
        y = rng.normal(d.mean(), d.std(ddof=1), n)
        o = np.argsort(y)
        sh.append(y[o[:10]].sum() / y.sum() if y.sum() != 0 else np.nan)
        a10.append(y[o[10:]].mean())
        a20.append(y[o[20:]].mean() if n > 20 else np.nan)
    sh = np.asarray(sh); a10 = np.asarray(a10); a20 = np.asarray(a20)
    q = lambda v: [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))]
    return dict(share=float(np.nanmean(sh)), share_ci=q(sh), share_draws=sh,
                dt10=float(np.nanmean(a10)), dt10_ci=q(a10), dt10_draws=a10,
                dt20=float(np.nanmean(a20)), dt20_ci=q(a20), dt20_draws=a20)


def check(path=None, k=30.0):
    """The full discipline on one restraint arm: drop-top-10/20, per fold, FAIL18.

    A -0.025 A mean is small enough that a handful of targets could carry it -- the
    `leg_torsion` lead this sprint looked like the best native-free selector in the project
    and was two targets of nine (E7-REFUTED).  This is the check that decides it.
    """
    import json
    import os
    p = path or os.path.join(E.RESULTS, "ener_avgrefine.json")
    with open(p) as fh:
        d = json.load(fh)
    rows = [r for r in d["per_target"] if f"amber_k{k}" in r and "rmsd" in r[f"amber_k{k}"]]
    A = np.array([r["projected_rmsd"] for r in rows])
    V = np.array([r[f"amber_k{k}"]["rmsd"] for r in rows])
    B = np.array([r["raw_avg_rmsd"] for r in rows])
    folds = [r["fold"] for r in rows]
    f18 = set(I.FAIL18)
    print(f"=== k={k:g} on {len(rows)} targets: A {A.mean():.3f}  B {B.mean():.3f}  "
          f"B+AMBER {V.mean():.3f} ===")
    p1 = I.paired(V, A, folds=folds)
    print(f"\nvs arm A (the incumbent projection)")
    print(f"  mean diff      {p1['mean_diff']:+.4f}  CI [{p1['ci95'][0]:+.4f},"
          f"{p1['ci95'][1]:+.4f}]  {p1['n_better']}W/{p1['n_worse']}L  "
          f"median {p1['median_diff']:+.4f}")
    d0 = V - A
    order = np.argsort(d0)
    dt10 = p1["drop_top10_mean_diff"]; dt20 = p1["drop_top20_mean_diff"]
    share = p1["top10_share"]
    print(f"\n  CONCENTRATION -- the three statistics that decide it, together:")
    print(f"    drop-top-10  {dt10:+.4f}      drop-top-20  {dt20:+.4f}      "
          f"top-10 share  {share:.3f}")
    # A raw drop-top threshold is NOT a valid concentration test.  When the effect is
    # small relative to the per-target spread, removing the 10 most favourable targets
    # mechanically eats most of the mean even if the effect is perfectly uniform.  So the
    # observed statistic must be compared with its own null: a uniform effect of the SAME
    # mean and the SAME per-target sd.  (This replaces an invented 0.5 threshold that
    # produced a false FAILS verdict -- see E14e.)
    nl = concentration_null(d0, seed=0)
    print(f"    NULL (a UNIFORM effect of the same mean and sd would give):")
    print(f"      drop-top-10  {nl['dt10']:+.4f} [{nl['dt10_ci'][0]:+.4f},"
          f"{nl['dt10_ci'][1]:+.4f}]   drop-top-20  {nl['dt20']:+.4f} "
          f"[{nl['dt20_ci'][0]:+.4f},{nl['dt20_ci'][1]:+.4f}]   "
          f"share {nl['share']:.3f} [{nl['share_ci'][0]:.3f},{nl['share_ci'][1]:.3f}]")
    pct20 = float((nl["dt20_draws"] < dt20).mean())
    pcts = float((nl["share_draws"] < share).mean())
    print(f"      observed sits at percentile {pct20:.2f} (drop-top-20) and "
          f"{pcts:.2f} (share) of that null")
    conc = (pct20 > 0.90) or (pcts > 0.90)
    print(f"    VERDICT: {'CONCENTRATED' if conc else 'NOT concentrated'} -- the "
          f"concentration statistics are "
          f"{'ABOVE' if conc else 'INDISTINGUISHABLE FROM'} what a uniform effect of this "
          f"size and noise produces.")
    print(f"    mean/sd = {d0.mean()/d0.std(ddof=1):+.3f}; when this ratio is small the "
          f"drop-top test has little power and cannot settle concentration either way.")
    print(f"    (drop-top-10 removes {10/len(rows)*100:.0f}% of this sample)")
    print(f"  decay ladder:  " + "  ".join(
        f"top{kk}:{d0[order[kk:]].mean():+.4f}" for kk in (1, 2, 3, 5)))
    print(f"\nper fold:")
    for f in sorted(set(folds)):
        m = np.array([x == f for x in folds])
        print(f"  fold {f}  n={m.sum():3d}  A {A[m].mean():.3f}  B+AMBER {V[m].mean():.3f}"
              f"  diff {d0[m].mean():+.4f}  "
              f"{int((d0[m]<0).sum())}W/{int((d0[m]>0).sum())}L")
    print(f"\nFAIL18 vs the rest:")
    for nm, m in (("FAIL18", np.array([r["pdb"] in f18 for r in rows])),
                  ("other", np.array([r["pdb"] not in f18 for r in rows]))):
        if m.sum():
            print(f"  {nm:8s} n={m.sum():3d}  A {A[m].mean():.3f}  "
                  f"B+AMBER {V[m].mean():.3f}  diff {d0[m].mean():+.4f}  "
                  f"{int((d0[m]<0).sum())}W/{int((d0[m]>0).sum())}L")
    g = np.array([r[f"amber_k{k}"]["geom"]["geom_rms_rel_dev"] for r in rows])
    print(f"\ngeometry: deviation {g.mean():.4f} (max {g.max():.4f})   "
          f"clashed {np.mean([r[f'amber_k{k}']['geom']['n_clash_2A'] > 0 for r in rows]):.3f}"
          f"   rama {np.mean([r[f'amber_k{k}']['rama_ok'] for r in rows]):.3f}")
    print(f"the 5 largest gains: " + ", ".join(
        f"{rows[i]['pdb']} {d0[i]:+.3f}" for i in order[:5]))
    print(f"the 5 largest losses: " + ", ".join(
        f"{rows[i]['pdb']} {d0[i]:+.3f}" for i in order[-5:]))
    return p1


def report_only(ks=(10.0, 2.0, 0.0)):
    """Re-report from the checkpointed JSON, so a partial run is still analysable."""
    import json
    import os
    with open(os.path.join(E.RESULTS, "ener_avgrefine.json")) as fh:
        d = json.load(fh)
    rows = d["per_target"]
    print(f"loaded {len(rows)} targets (complete={d.get('complete')})")
    report(rows, tuple(d.get("ks", ks)))
    return rows


if __name__ == "__main__":
    a = sys.argv[1:]
    if "--check" in a:
        kk = float(a[a.index("--check") + 1]) if len(a) > a.index("--check") + 1             and not a[a.index("--check") + 1].startswith("-") else 30.0
        pth = a[a.index("--path") + 1] if "--path" in a else None
        check(pth, kk)
    elif "--report-only" in a:
        report_only()
    else:
        pil = int(a[a.index("--pilot") + 1]) if "--pilot" in a else None
        kk = tuple(float(x) for x in a[a.index("--k") + 1].split(",")) \
            if "--k" in a else (10.0, 2.0, 0.0)
        main(pil, kk)
