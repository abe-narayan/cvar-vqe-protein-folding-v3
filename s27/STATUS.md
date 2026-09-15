
# S28 STATUS

## Coordinator
- 2026-09-14 19:00 running: governor up, cap 4; spawning lanes A (amplitude readout), B (hopping Hamiltonian), C (FAIL18 detector + readouts), D (Adversary).
- 2026-09-14 19:05 running: four lanes spawned 19:02 (A, B, C, D). Box: user baseline ~10 GB (chrome 4.5, claude 2.1, Code 1.9); RAM 88 to 91% before any lane job; jobrun will serialise jobs needing more than ~1 GB. Lanes: keep est-ram honest, checkpoint every 10 min, expect waits.
- next: hourly ledger commits; reassign lanes as they finish; final S28 report in s27/.

## A

## B

## C
- 2026-09-14 19:12 running: contract, brief, S27 report/prereg, router record (S22 L7, S23 L6/L7, S26 L110/L115) read; PREREG_S28_C.md written (no endpoint read). Next: s28_C_fail18.py features + nested logistic + permutation null; s28_C_readout.py; tests; 1-target probe under jobrun.

## D
