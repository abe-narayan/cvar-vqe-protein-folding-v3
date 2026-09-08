"""s18/phys_space.py -- PHASE 10 (SECONDARY): CHEAPER SPACING RESTORATION.

Pre-registration: `s18/PREREG_phys.md` section 4.  This phase is explicitly secondary and is not
allowed to crowd out Phases 6 and 8; it is therefore split so that the CHEAP half can falsify it
without spending a single AMBER minimisation.

THE DEFECT, from this workstream's own Sprint-17 measurement (ledger L28).  The coordinate
average contracts Ca-Ca spacing to **2.949 A against an ideal trans value of 3.80 A** -- 22.4% --
and

    Spearman(input Ca-Ca spacing, cis fraction at `cafix`)   = -0.949
    Spearman(input Ca-Ca spacing, convergence at `cafix`)    = +0.807

so the peptide-bond defect and the 38% gate-failure rate are both almost deterministic functions
of that one number.  The incumbent fix -- the ideal-geometry projection -- restores the spacing
and costs **+0.164 A [+0.131, +0.199], 15W/111L** (L30).  The question is whether a cheaper
native-free correction recovers part of that tax.

ARMS, all native-free, all applied to the SAME all-atom coordinate average:

    none            the 3.048 A average                                  the reference
    scale_uniform   ONE parameter: rigid dilation about the centroid to mean spacing 3.80
    bond_norm       LOCAL: every Ca-Ca bond direction kept, every length set to 3.80
    segment         SEGMENT-WISE: bond lengths rescaled by a windowed local target
    minproj         MINIMAL DISPLACEMENT: min ||X - X0||^2 s.t. |x_{i+1}-x_i| = 3.80 (SHAKE)
    torsion_rebuild CONSTRAINED RECONSTRUCTION preserving the AVERAGED TORSIONS exactly
    proj            the incumbent ideal-geometry projection               the +0.164 A reference
    proj_fit        its unpenalised (lambda = 0) twin
    CONTROLS
    helix           ZERO-INFORMATION: the constant ideal alpha-helix (spacing 3.80 by construction)
    rand_move       MATCHED-RANDOM displacement of each arm's own realised magnitude

PASS 1 (cheap, n = 126, no AMBER): realised spacing and Ca-RMSD cost.  If no arm restores the
spacing at less Ca cost than the projection, **the branch closes here** and no AMBER time is
spent.  PASS 2 (AMBER k = 30) runs only on the arms Pass 1 leaves alive.

WHY A VALIDITY NUMBER ALONE WOULD BE MEANINGLESS HERE, and it is the whole point of the two
controls.  A constant alpha-helix has spacing 3.80, Ramachandran 1.000 and zero clashes **by
construction** and is 4.07 A from the native.  Restoring the spacing is trivial; restoring it
without paying more than 0.164 A of accuracy is the entire question.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from s12 import instrument as I                                   # noqa: E402
from s14.avgspace import top75_windows                            # noqa: E402
from s14.ener_avgrefine import torsions_of                        # noqa: E402
from s15.phys_repl import averaged_backbone_from, project_both    # noqa: E402
from s15 import seed as SD                                        # noqa: E402
from s16 import energy_lib as EL                                  # noqa: E402
from s17 import phys_lib as P17                                   # noqa: E402
from s18 import phys_lib as PL                                    # noqa: E402

ATOMS = ("N", "CA", "C", "O", "CB")
IDEAL_CA = 3.80                  # trans Ca-Ca, the value L28's 22.4% is measured against
SEG_W = 3                        # window for the segment-wise arm
SHAKE_SWEEPS = 200
K_RESTRAINT = 30.0


def spacing(CA):
    return float(np.linalg.norm(np.diff(np.asarray(CA, float), axis=0), axis=1).mean())


# ------------------------------------------------------------------ the corrections
def scale_uniform(CA):
    """ONE parameter, native-free: dilate about the centroid until mean spacing is 3.80."""
    X = np.asarray(CA, float)
    c = X.mean(0)
    f = IDEAL_CA / max(spacing(X), 1e-9)
    return c + (X - c) * f


def bond_norm(CA, target=None):
    """Keep every Ca-Ca bond DIRECTION, set every length to its target (default 3.80)."""
    X = np.asarray(CA, float)
    b = np.diff(X, axis=0)
    L = np.linalg.norm(b, axis=1, keepdims=True)
    tg = np.full((len(b), 1), IDEAL_CA) if target is None else np.asarray(target, float).reshape(-1, 1)
    out = np.empty_like(X)
    out[0] = X[0]
    out[1:] = X[0] + np.cumsum(b / np.maximum(L, 1e-9) * tg, axis=0)
    #: superpose back onto the input so the correction is a SHAPE change, not a rigid move
    return I.superpose_batch(out[None], X)[0]


def segment(CA, w=SEG_W):
    """Segment-wise: each bond is rescaled toward 3.80 by the ratio its own WINDOW needs.

    Between `scale_uniform` (one global factor) and `bond_norm` (every bond forced exactly):
    a bond that is locally near-ideal is left nearly alone.
    """
    X = np.asarray(CA, float)
    L = np.linalg.norm(np.diff(X, axis=0), axis=1)
    k = np.ones(w) / w
    Ls = np.convolve(np.pad(L, (w // 2, w - 1 - w // 2), mode="edge"), k, mode="valid")
    return bond_norm(X, target=L * (IDEAL_CA / np.maximum(Ls, 1e-9)))


def minproj(CA, sweeps=SHAKE_SWEEPS):
    """MINIMAL DISPLACEMENT onto the spacing constraint set, by SHAKE-style projection.

    min ||X - X0||^2 subject to |x_{i+1} - x_i| = 3.80.  Each sweep corrects every bond by
    moving BOTH endpoints equally and oppositely, which is the steepest-descent projection of
    the single constraint; sweeping to convergence is Dykstra's algorithm on a chain.
    """
    X = np.asarray(CA, float).copy()
    for _ in range(sweeps):
        mx = 0.0
        for i in range(len(X) - 1):
            b = X[i + 1] - X[i]
            L = float(np.linalg.norm(b))
            if L < 1e-9:
                continue
            corr = 0.5 * (IDEAL_CA - L) / L
            X[i] -= corr * b
            X[i + 1] += corr * b
            mx = max(mx, abs(IDEAL_CA - L))
        if mx < 1e-6:
            break
    return I.superpose_batch(X[None], np.asarray(CA, float))[0]


def lift_ca(bb, CA_new):
    """Lift a Ca-only correction to all atoms by translating each residue rigidly.

    Every intra-residue bond, angle and chirality is preserved EXACTLY; only the inter-residue
    spacing changes.  This is what makes a Ca-only correction usable by AMBER without paying for
    the projection it is supposed to replace.
    """
    d = np.asarray(CA_new, float) - np.asarray(bb["CA"], float)
    return {a: np.asarray(bb[a], float) + d for a in bb}


def torsion_rebuild(avg, seq):
    """CONSTRAINED RECONSTRUCTION preserving the AVERAGED TORSIONS exactly.

    Read (phi, psi) off the averaged all-atom backbone and rebuild it with ideal bond geometry.
    Unlike the projection this fits nothing and has no Ramachandran penalty: it is exactly the
    averaged torsions on an ideal chain.
    """
    from core import geometry as geo
    ph, ps = torsions_of({a: np.asarray(avg[a], float) for a in ATOMS if a in avg})
    bb = geo.build_backbone(np.asarray(ph, float), np.asarray(ps, float))
    ca = I.superpose_batch(np.asarray(bb["CA"], float)[None], np.asarray(avg["CA"], float))[0]
    return ca, bb, np.asarray(ph, float), np.asarray(ps, float)


# ------------------------------------------------------------------ one target
def run_target(t, want_amber=False):
    pdb, seq, fold, n = t["pdb"], t["seq"], int(t["fold"]), int(t["n"])
    W, PHI, PSI, u = top75_windows(pdb)
    nat = np.asarray(u["nat_ca"], float)
    avg, C_ca, _dev = averaged_backbone_from(np.asarray(W, float), PHI, PSI)
    X0 = np.asarray(avg["CA"], float)
    rng = SD.stable_rng(pdb, "s18space")

    rec = {"pdb": pdb, "n": n, "fold": fold, "arms": {}}

    def add(nm, ca, bb=None):
        ca = np.asarray(ca, float)
        r = {"rmsd": float(I.ca_rmsd(ca, nat)), "spacing": spacing(ca),
             "disp": float(I.ca_rmsd(ca, X0))}
        if bb is not None:
            r["valid"] = EL.panel({a: np.asarray(bb[a], float) for a in ATOMS if a in bb}, seq)
        rec["arms"][nm] = r
        return r

    add("none", X0, avg)
    add("scale_uniform", scale_uniform(X0))
    add("bond_norm", bond_norm(X0))
    add("segment", segment(X0))
    add("minproj", minproj(X0))
    ca_tr, bb_tr, ph_tr, ps_tr = torsion_rebuild(avg, seq)
    add("torsion_rebuild", ca_tr, bb_tr)

    from core import geometry as geo
    pr = project_both(C_ca, seq, fold)
    add("proj", pr["ca"], geo.build_backbone(pr["phi"], pr["psi"]))
    add("proj_fit", pr["fit_ca"], geo.build_backbone(pr["fit_phi"], pr["fit_psi"]))

    #: CONTROL -- ZERO INFORMATION.  Spacing 3.80 and Ramachandran 1.000 BY CONSTRUCTION.
    hb = P17.helix_backbone(n)
    hca = I.superpose_batch(np.asarray(hb["CA"], float)[None], X0)[0]
    add("helix", hca, hb)

    #: CONTROL -- MATCHED RANDOM, one per arm at that arm's OWN realised displacement.
    for nm in ("scale_uniform", "bond_norm", "segment", "minproj", "torsion_rebuild", "proj"):
        z = P17.matched_random_disp(X0, rec["arms"][nm]["disp"], rng)
        rec["arms"]["rand@" + nm] = {"rmsd": float(I.ca_rmsd(z, nat)),
                                     "spacing": spacing(z),
                                     "disp": float(I.ca_rmsd(z, X0))}

    if want_amber:
        from core import amber as am
        import torsion_lib2 as tl2
        tab = tl2.library_for(seq, 4, seq)
        rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
        from core import geometry as geo2
        #: PASS 2 runs on the arms PASS 1 leaves alive: the un-corrected average (the input whose
        #: contraction causes the defect), the minimal-displacement Ca correction, and the
        #: incumbent projection -- the +0.164 A reference the cheaper arm has to beat.
        #:
        #: THE CA-ONLY ARMS HAVE TO BE LIFTED TO ALL ATOMS, and how is a real choice, not a
        #: detail.  Re-projecting them would defeat the point (the projection is the expensive
        #: thing being replaced), so each residue's whole atom group is TRANSLATED rigidly by its
        #: own Ca displacement.  That preserves every intra-residue bond and angle exactly and
        #: changes only the inter-residue spacing, which is the quantity under test.
        ph_p, ps_p = np.asarray(pr["phi"], float), np.asarray(pr["psi"], float)
        for nm, bb in (("none", avg), ("torsion_rebuild", bb_tr),
                       ("minproj", lift_ca(avg, minproj(X0))),
                       ("proj", geo2.build_backbone(ph_p, ps_p))):
            coords = {a: np.asarray(bb[a], float) for a in ATOMS if a in bb}
            t0 = time.time()
            r = am.refine_coords(seq, rep, coords, k_restraint=K_RESTRAINT, steps=0,
                                 tolerance=1.0, threads=1, memo=False)
            ca = np.asarray(r["ca"], float)
            ob = {a: np.asarray(r["backbone"][a], float) for a in ATOMS if a in r["backbone"]}
            rec["arms"][nm]["amb30"] = {
                "rmsd": float(I.ca_rmsd(ca, nat)),
                "input_rmsd": float(I.ca_rmsd(coords["CA"], nat)),
                "energy": float(r["energy"]), "converged": bool(r["converged"]),
                "ca_disp": float(I.ca_rmsd(ca, coords["CA"])),
                "valid": EL.panel(ob, seq), "input_valid": EL.panel(coords, seq),
                "wall": round(time.time() - t0, 2)}
    return rec


def run(targets=None, out="space.json", want_amber=False, verbose=True):
    tg = targets if targets is not None else I.targets()
    cfg = {"IDEAL_CA": IDEAL_CA, "SEG_W": SEG_W, "SHAKE_SWEEPS": SHAKE_SWEEPS,
           "want_amber": bool(want_amber), "K_RESTRAINT": K_RESTRAINT}
    path = os.path.join(PL.RESULTS, out)
    rows = []
    if os.path.exists(path):
        try:
            prev = json.load(open(path))
            if prev.get("cfg_hash") == PL.cfg_hash(cfg):
                rows = prev.get("rows", [])
        except Exception:
            rows = []
    done = {r["pdb"] for r in rows}
    t0 = time.time()
    for t in tg:
        if t["pdb"] in done:
            continue
        PL.mem_ok(1.2)
        rows.append(run_target(t, want_amber=want_amber))
        if verbose and len(rows) % 10 == 0:
            print(f"  {len(rows)}/{len(tg)} ({time.time()-t0:.0f}s)", flush=True)
        PL.write(out, {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg)},
                 n_expected=len(tg))
    PL.write(out, {"rows": rows, "config": cfg, "cfg_hash": PL.cfg_hash(cfg)}, n_expected=len(tg))
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--amber", action="store_true")
    ap.add_argument("--out", default="space.json")
    a = ap.parse_args()
    tg = I.targets()
    if a.limit:
        tg = tg[:a.limit]
    run(tg, out=a.out, want_amber=a.amber)
