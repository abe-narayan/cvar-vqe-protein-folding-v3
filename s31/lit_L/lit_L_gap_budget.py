"""Does the VQE-vs-closed-form gap close with more budget, or is it an EXPRESSIVITY floor?

The ansatz has n*layers parameters for a (2**n - 1)-dimensional simplex: 21 params for 127
dimensions at n=7,layers=3.  If the gap plateaus above zero as iters/restarts/layers grow,
the circuit cannot EXPRESS the closed-form optimum p*, and no optimiser fixes it.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, r"C:\Users\abena\Protein-Folding-Algorithm")
from core.quantum import cvar_exact, run_cvar_vqe

sys.path.insert(0, os.path.dirname(__file__))
from lit_L_dequant_check import closed_form, F_of_p


def main():
    n = 7
    D = 1 << n
    r = np.arange(D, dtype=float)
    E0 = (r - r.mean()) / r.std()
    rng = np.random.default_rng(11)
    rng.shuffle(E0)

    print(f"{'alpha':>5} {'T':>5} {'layers':>6} {'iters':>6} {'rst':>4} |"
          f" {'F_vqe':>11} {'F_closed':>11} {'gap':>10} {'TV':>6} {'H_vqe':>7} {'sec':>6}")
    print("-" * 96)
    for alpha, T in ((0.1, 0.1), (1.0, 0.1)):
        p_cf, F_dual, s = closed_form(E0, alpha, T)
        F_cf, _, H_cf = F_of_p(E0, p_cf, alpha, T)
        for layers in (3, 6, 12):
            for iters, rst in ((80, 1), (400, 1), (2000, 1), (400, 8)):
                t0 = time.time()
                p_q, cv, H, _ = run_cvar_vqe(E0, alpha, T, n=n, layers=layers,
                                             iters=iters, restarts=rst, seed=0)
                dt = time.time() - t0
                F_q = float(cv - T * H)
                tv = 0.5 * float(np.abs(p_q - p_cf).sum())
                print(f"{alpha:5.2f} {T:5.2f} {layers:6d} {iters:6d} {rst:4d} |"
                      f" {F_q:11.6f} {F_cf:11.6f} {F_q - F_cf:+10.6f} {tv:6.3f}"
                      f" {H:7.3f} {dt:6.1f}")
        print(f"{'':>5} {'':>5} {'closed form:':>44} H_cf = {H_cf:.3f}, s* = {s:.4f}")
        print()


if __name__ == "__main__":
    main()
