#!/usr/bin/env python
"""s27/s28_B_split.py -- the THREE-WAY SPLIT's middle leg (S28-L2 caveat (b)).

For every (target, graph, J) of the S28-B run: the objective F = CVaR_alpha(E; p) - T H(p)
- J <psi|A|psi> (alpha 0.18, T 0.5, the objective the circuit optimises) evaluated at

    F_vqe        the VQE optimum (read from `s28_B_rows.jsonl`, seeds 0 and 1)
    F_gs         the exact ground state of H = diag(E) - J A (a state the circuit could in
                 principle represent; it is the optimum of <psi|H|psi>, i.e. of F at alpha = 1,
                 T = 0, NOT of the circuit's objective)
    F_untrained  the best of 16 untrained draws theta ~ N(0, 0.6^2) (draw 1 is the VQE's own
                 initial theta for that seed, so F_untrained <= F(theta_0))

and the realised hopping value <psi|A|psi> at each, against its exact maximum 1.0 (unit
spectral norm; the Perron vector attains it). Native-free; no RMSD. Writes `s28_B_split.json`.
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
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from core import quantum as Q              # noqa: E402
from s22 import qcand_lib as QC            # noqa: E402
from s24 import stats_lib as ST            # noqa: E402
from s27 import run_pool as RP             # noqa: E402
from s27 import s28_B_hop as B             # noqa: E402

OUT = os.path.join(B.RESULTS, "s28_B_split.json")
ROWS = os.path.join(B.RESULTS, "s28_B_split_rows.jsonl")
N_UNTRAINED = 16


def F_of_state(psi: np.ndarray, E: np.ndarray, A: np.ndarray, J: float) -> dict:
    """The circuit's objective at an arbitrary real unit state (free_energy's arithmetic)."""
    p = psi ** 2
    p = p / p.sum()
    v, _, _ = Q.cvar_exact(E, p, B.ALPHA)
    lp = np.log(np.maximum(p, 1e-15))
    Hn = float(-(p * lp).sum())
    hv = B.hop_value(psi, A)
    return dict(F=float(v - B.TEMP * Hn - J * hv), cvar=float(v), H=Hn, hop=float(hv),
                sign_coh=B.sign_coherence(psi))


def main():
    rows = B.load_rows()
    by = {}
    for r in rows:
        by.setdefault(r["pdb"], {})[B.arm_of(r)] = r
    done = B._done_pdbs(ROWS)
    t0 = time.time()
    for pdb in sorted(by):
        if pdb in done:
            continue
        cand, ch, _ = RP.channels_for(pdb)
        E = RP.zr(ch["DIS"])
        enc = QC.Encoding(E)
        D, G = B.build_graphs(cand.W, pdb)
        circ = Q.StatevectorCircuit(enc.n_qubits, B.LAYERS)
        P = circ.n_params()
        out = []
        for gname in B.GRAPHS:
            A_pad = B.pad_graph(G[gname]["A"], enc.dim)
            for J in B.J_GRID:
                if J == 0.0 and gname != "REAL":
                    continue
                A_use = np.zeros_like(A_pad) if J == 0.0 else A_pad
                gtag = "NONE" if J == 0.0 else gname
                gs = B.ground_state(enc.E, A_use, J)
                f_gs = F_of_state(gs["psi"], enc.E, A_use, J)
                rec = dict(pdb=pdb, graph=gtag, J=float(J), e0=gs["e0"], F_gs=f_gs["F"],
                           hop_gs=f_gs["hop"], cvar_gs=f_gs["cvar"], H_gs=f_gs["H"])
                for s in B.SEEDS:
                    rng = np.random.default_rng(s)
                    fs = []
                    for _ in range(N_UNTRAINED):
                        th = rng.normal(0.0, 0.6, P)
                        fs.append(F_of_state(circ.state(th), enc.E, A_use, J))
                    fb = min(fs, key=lambda d: d["F"])
                    vq = by[pdb][f"vqe|s{s}|{gtag}|J{J:g}"]
                    rec[f"F_vqe_s{s}"] = float(vq["F"])
                    rec[f"hop_vqe_s{s}"] = float(vq.get("hop", 0.0))
                    rec[f"coh_vqe_s{s}"] = float(vq.get("sign_coh", 0.0))
                    rec[f"F_untrained_best16_s{s}"] = float(fb["F"])
                    rec[f"F_untrained_theta0_s{s}"] = float(fs[0]["F"])
                    rec[f"hop_untrained_best16_s{s}"] = float(fb["hop"])
                out.append(rec)
        with open(ROWS, "a", encoding="utf-8") as fh:
            for r in out:
                fh.write(json.dumps(r) + "\n")
        print(f"  {pdb} ({(time.time()-t0)/60:.1f} min)", flush=True)
    srows = B.load_rows(ROWS)
    summ = {}
    for gname in ("NONE",) + tuple(B.GRAPHS):
        for J in B.J_GRID:
            rs = [r for r in srows if r["graph"] == gname and r["J"] == J]
            if not rs:
                continue
            d = dict(n=len(rs))
            for k in ("F_gs", "hop_gs", "e0"):
                d[k + "_mean"] = float(np.mean([r[k] for r in rs]))
            for s in B.SEEDS:
                fv = np.array([r[f"F_vqe_s{s}"] for r in rs]); fg = np.array([r["F_gs"] for r in rs])
                fu = np.array([r[f"F_untrained_best16_s{s}"] for r in rs])
                d[f"F_vqe_s{s}_mean"] = float(fv.mean())
                d[f"F_untrained_best16_s{s}_mean"] = float(fu.mean())
                d[f"vqe_below_gs_frac_s{s}"] = float(np.mean(fv < fg))
                d[f"vqe_below_untrained_frac_s{s}"] = float(np.mean(fv < fu))
                d[f"gap_closed_s{s}"] = float(np.mean((fu - fv) / np.maximum(fu - np.minimum(fg, fv), 1e-12)))
                d[f"hop_vqe_s{s}_mean"] = float(np.mean([r[f"hop_vqe_s{s}"] for r in rs]))
                d[f"coh_vqe_s{s}_mean"] = float(np.mean([r[f"coh_vqe_s{s}"] for r in rs]))
            summ[f"{gname}|J{J:g}"] = d
    ST.save_atomic(OUT, dict(rows=srows, summary=summ, n_untrained=N_UNTRAINED), module_file=__file__)
    print(f"\n{'cell':12s} {'F_vqe s0':>9} {'F_vqe s1':>9} {'F_gs':>9} {'F_untr16 s0':>11} {'vqe<gs s0':>9} {'vqe<untr s0':>11} {'hop_vqe s0':>10} {'hop_gs':>7} {'coh s0':>6}")
    for k, d in summ.items():
        print(f"{k:12s} {d['F_vqe_s0_mean']:9.4f} {d['F_vqe_s1_mean']:9.4f} {d['F_gs_mean']:9.4f} {d['F_untrained_best16_s0_mean']:11.4f} "
              f"{d['vqe_below_gs_frac_s0']:9.3f} {d['vqe_below_untrained_frac_s0']:11.3f} {d['hop_vqe_s0_mean']:10.3f} {d['hop_gs_mean']:7.3f} {d['coh_vqe_s0_mean']:6.3f}")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
