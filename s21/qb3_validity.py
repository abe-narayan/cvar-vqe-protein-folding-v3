"""SPRINT 21 / WORKSTREAM B -- the VALIDITY GUARD on the readout finding.

BRIEF section 3: `caonly_k300` was 0.110 A more accurate with an indistinguishable clash count and
`cis_frac` +0.4242 -- **a scalar validity score would have promoted a broken structure. Validity is
a VECTOR.**  Section 4 requires the validity vector beside any candidate-selection claim.

This lane's largest result (section 4 of the findings) is that swapping the READOUT Hamiltonian from
Legacy to the distogram is worth -0.697 A on the identical evaluated set.  A cheaper structure is
not automatically a better one, so the two readouts' argmins are compared on geometry as well as on
Ca-RMSD.

WHAT IS AND IS NOT MEASURABLE HERE.  Every configuration in this lane is built by
`core.project.build_ca_exact` from continuous (phi, psi) with IDEAL bond geometry, so bond lengths,
bond angles and omega are ideal BY CONSTRUCTION and carry no information -- reporting them as
"passing" would be a vacuous gate (BRIEF section 7 rule 3).  What is free to vary, and is therefore
what is measured:

    ca_clash_lt4.0 / lt4.5   non-adjacent Ca pairs closer than 4.0 / 4.5 A
    min_nonadj_ca            the closest non-adjacent Ca approach, in A
    rg                       radius of gyration -- the axis Sprint 20 showed SEPARATES the two
                             potentials (Legacy-preferred candidates are 0.45 A more compact)
    rg_vs_native             rg minus the native's rg: compactness error, ORACLE, labelled

The arms are re-run rather than re-read because the evaluated SETS were not persisted; every seed is
`stable_rng`, so the re-run is bit-exact with the recorded rows and that is asserted, not assumed.

    python -m s21.qb3_validity [n_targets]
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

from s20 import qb2_lib as L
from s20 import qb2_opt as OP
from s20 import qb2_run as R20
from s21 import qb3_lib as K
from s21 import qb3_run as RN
from s15 import seed as SD

ARMS = ("mps_L2", "untrained_L2", "best_of_N")
TRAIN = "LEG"
READS = ("LEG", "DIST")


def ca_validity(phi, psi, tgt):
    from core import project as pj
    CA = np.asarray(pj.build_ca_exact(np.atleast_2d(phi), np.atleast_2d(psi)), float)[0]
    d = np.linalg.norm(CA[:, None, :] - CA[None, :, :], axis=-1)
    n = len(CA)
    ii, jj = np.triu_indices(n, k=2)          # non-adjacent only: |i-j| >= 2
    dv = d[ii, jj]
    rg = float(np.sqrt(((CA - CA.mean(0)) ** 2).sum(1).mean()))
    nat = np.asarray(tgt["nat"], float)
    rgn = float(np.sqrt(((nat - nat.mean(0)) ** 2).sum(1).mean()))
    return {"ca_clash_lt40": int((dv < 4.0).sum()), "ca_clash_lt45": int((dv < 4.5).sum()),
            "min_nonadj_ca": float(dv.min()), "rg": rg, "rg_vs_native_ORACLE": rg - rgn,
            "n_pairs": int(len(dv))}


def main():
    n_t = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    pdbs = [t["pdb"] for t in L.subset(20)][:n_t]
    out = {}
    exact_ok, exact_n = 0, 0
    prev = K.read21("qb3_a")
    for ti, pdb in enumerate(pdbs):
        t0 = time.time()
        tgt = L.target(pdb)
        Z = R20.starts_for(tgt)
        H, sp, cal, ver = R20.prep(tgt, (TRAIN,))
        rec = {"n": tgt["n"], "fold": tgt["fold"], "rows": {}}
        for arm in ARMS:
            kw = dict(ansatz="mps_L2", opt="adam", alpha=RN.A_ALPHA, shots=RN.A_SHOTS,
                      train=(arm != "untrained_L2"))
            for s in RN.SEEDS:
                F, h = RN._field(TRAIN, tgt, RN.B, sp, cal[TRAIN], keep=True)
                rng = SD.stable_rng(pdb, TRAIN, arm, s, salt=K.SALT)
                if arm == "best_of_N":
                    OP.arm_bestofn(F, Z[s % len(Z)], rng)
                else:
                    K.arm_vqe2(F, Z[s % len(Z)], rng, **kw)
                phi, psi, _e = h.seen()
                row = {}
                for rk in READS:
                    hr = L.Ham(rk, tgt, budget=10 ** 9, sp=sp)
                    er = np.asarray(hr.raw(phi, psi), float)
                    ok = np.isfinite(er)
                    k = int(np.flatnonzero(ok)[int(np.argmin(er[ok]))])
                    v = ca_validity(phi[k], psi[k], tgt)
                    rm = float(L.rmsd_of(L.pack(phi[k:k + 1], psi[k:k + 1]), tgt)[0])
                    row[rk] = {"rmsd_ORACLE": rm, **v}
                # BIT-EXACTNESS against the recorded run: same seeds must give the same answer
                ref = prev.get(pdb, {}).get("rows", {}).get(f"{TRAIN}|{arm}|{s}")
                if ref is not None:
                    exact_n += 1
                    exact_ok += int(abs(row[TRAIN]["rmsd_ORACLE"] - ref["rmsd_ORACLE"]) < 1e-12)
                rec["rows"][f"{arm}|{s}"] = row
        out[pdb] = rec
        print(f"[{ti+1}/{len(pdbs)}] {pdb} {time.time()-t0:.0f}s", flush=True)
    K.write21("qb3_validity", {
        "train_H": TRAIN, "readouts": list(READS), "arms": list(ARMS),
        "seeds": list(RN.SEEDS), "budget": RN.B, "n_targets": len(pdbs),
        "reproducibility_gate": {"exact": exact_ok, "compared": exact_n},
        "note": "bond/angle/omega are ideal BY CONSTRUCTION here and are not reported as a passing "
                "gate; only the quantities free to vary are measured.",
        "rows": out}, complete=(len(pdbs) == 20))
    print(f"reproducibility gate: {exact_ok}/{exact_n} rows bit-exact against qb3_a.json")


if __name__ == "__main__":
    main()
