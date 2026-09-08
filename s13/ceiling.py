"""SPRINT 13, EXPERIMENT 0 -- the representation ceiling of the discrete torsion space.

This is the gating number for the whole sprint.  Sprint 12's oracle diagnostic showed that
CONTINUOUS backbone torsions at sigma = 12 degrees build a 1.486 A structure.  The proposed
architecture optimises over a DISCRETE, sequence-conditioned torsion library instead --
`torsion_lib2.library_for(seq, k, exclude_seq=seq)` gives k states per residue, which is
`n_res * log2(k)` qubits.  Those are not the same object, and the discretisation costs
something.  If k=8 can only express 2.5 A structures then no optimiser, quantum or
classical, can do better than 2.5 A and the architecture is capped before it is built.

Three quantities per target, all ORACLE/DIAGNOSTIC (they read native torsions):

    snap      assign each residue independently to its nearest library state, build, measure.
              Greedy and chain-blind: it ignores that a small early torsion error is
              amplified by every downstream residue.  An upper bound on the ceiling.
    descent   coordinate descent over residues on the TRUE CA-RMSD, from the snap start.
              Chain-aware, so it recovers what greedy loses.  A tighter upper bound.
    best      the better of the two.

The ceiling is a property of the LIBRARY, and the library is built with the target's own
sequence excluded, so this is not leakage in the retrieval sense -- but every arm here reads
the native to score, so all three are diagnostics and none is a predictive result.

    python -m s13.ceiling [k ...]
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from s12 import instrument as I            # noqa: E402
import peptide_db as pdb                   # noqa: E402
import torsion_lib2 as tl2                 # noqa: E402

RESULTS = os.path.join(ROOT, "s13", "results")
os.makedirs(RESULTS, exist_ok=True)
KS = (4, 8, 16, 32)


def _wrap(a):
    return (np.asarray(a, float) + np.pi) % (2 * np.pi) - np.pi


def rep_tables(seq, k):
    """(n, k) phi and psi tables for the sequence-conditioned library, target held out."""
    tab = tl2.library_for(seq, k, seq)
    rep = tl2.PerResidueTorsion(seq, tab, chi_bits=False)
    return np.asarray(rep._phi, float), np.asarray(rep._psi, float), rep


def build(PHI, PSI, states):
    """CA trace for one per-residue state assignment."""
    n = len(states)
    phi = PHI[np.arange(n), states]
    psi = PSI[np.arange(n), states]
    return I.build_ca(phi, psi), phi, psi


def snap(PHI, PSI, phi0, psi0):
    """Nearest library state per residue in circular torsion distance. Chain-blind."""
    d = np.abs(_wrap(PHI - phi0[:, None])) + np.abs(_wrap(PSI - psi0[:, None]))
    return np.argmin(d, axis=1)


def descent(PHI, PSI, nat, states, max_sweeps=12):
    """Coordinate descent on the TRUE CA-RMSD. Chain-aware; the honest ceiling probe."""
    n, k = PHI.shape
    s = states.copy()
    best = I.ca_rmsd(build(PHI, PSI, s)[0], nat)
    for _ in range(max_sweeps):
        moved = False
        for i in range(n):
            cur = s[i]
            cand = np.repeat(s[None, :], k, axis=0)
            cand[:, i] = np.arange(k)
            W = np.stack([build(PHI, PSI, c)[0] for c in cand])
            r = I.kabsch_rmsd_batch(W, nat)
            b = int(np.argmin(r))
            if r[b] < best - 1e-9:
                best = float(r[b]); s[i] = b; moved = True
            else:
                s[i] = cur
        if not moved:
            break
    return s, best


def run(ks=KS):
    tg = I.targets()
    rows = []
    t0 = time.time()
    for idx, t in enumerate(tg):
        u = I.load_univ(t["pdb"]); nat = u["nat_ca"]
        p = pdb.by_pdb(t["pdb"])
        phi0 = np.asarray(p.phi, float); psi0 = np.asarray(p.psi, float)   # ORACLE
        row = {"pdb": t["pdb"], "n": t["n"], "fold": t["fold"],
               "native_rebuild": I.ca_rmsd(I.build_ca(phi0, psi0), nat)}
        for k in ks:
            PHI, PSI, _ = rep_tables(t["seq"], k)
            s0 = snap(PHI, PSI, phi0, psi0)
            r_snap = I.ca_rmsd(build(PHI, PSI, s0)[0], nat)
            s1, r_desc = descent(PHI, PSI, nat, s0)
            row[f"k{k}"] = {"snap": float(r_snap), "descent": float(r_desc),
                            "qubits": int(t["n"] * int(np.log2(k))),
                            "configs": float(k) ** t["n"],
                            "states": [int(x) for x in s1]}
        rows.append(row)
        if idx % 10 == 0 or idx == len(tg) - 1:
            json.dump({"what": "discrete torsion representation ceiling", "per_target": rows},
                      open(os.path.join(RESULTS, "ceiling.json"), "w"), indent=1)
            print(f"  {idx+1}/{len(tg)} {t['pdb']} n={t['n']} "
                  + "  ".join(f"k{k}:{row[f'k{k}']['snap']:.2f}->{row[f'k{k}']['descent']:.2f}"
                              for k in ks) + f"  [{time.time()-t0:.0f}s]", flush=True)
    json.dump({"what": "discrete torsion representation ceiling", "per_target": rows},
              open(os.path.join(RESULTS, "ceiling.json"), "w"), indent=1)
    report(rows, ks)
    return rows


def report(rows, ks=KS):
    f18 = set(I.FAIL18); isf = np.array([r["pdb"] in f18 for r in rows])
    nb = np.array([r["native_rebuild"] for r in rows], float)
    out = {"n": len(rows), "native_rebuild": I.summary(nb), "k": {}}
    print(f"\nideal-geometry rebuild of the NATIVE torsions: {nb.mean():.4f} "
          f"(the floor any torsion representation inherits)\n")
    hdr = (f"{'k':>3s} {'qubits':>7s} {'snap':>7s} {'descent':>8s} {'<2A':>6s} "
           f"{'<1.5A':>6s} {'FAIL18':>7s} {'other108':>8s}")
    print(hdr); print("-" * len(hdr))
    for k in ks:
        sn = np.array([r[f"k{k}"]["snap"] for r in rows], float)
        de = np.array([r[f"k{k}"]["descent"] for r in rows], float)
        qb = np.array([r[f"k{k}"]["qubits"] for r in rows], float)
        out["k"][k] = {"snap": I.summary(sn), "descent": I.summary(de),
                       "mean_qubits": float(qb.mean()), "max_qubits": float(qb.max()),
                       "FAIL18": float(de[isf].mean()), "other108": float(de[~isf].mean()),
                       "frac_under_1.5": float((de < 1.5).mean())}
        print(f"{k:3d} {qb.mean():7.1f} {sn.mean():7.3f} {de.mean():8.3f} "
              f"{(de<2.0).mean():6.2f} {(de<1.5).mean():6.2f} "
              f"{de[isf].mean():7.3f} {de[~isf].mean():8.3f}")
    with open(os.path.join(RESULTS, "ceiling_report.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {RESULTS}/ceiling_report.json")


if __name__ == "__main__":
    ks = tuple(int(x) for x in sys.argv[1:]) or KS
    run(ks)
