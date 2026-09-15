
# S28 STATUS

## Coordinator
- 2026-09-14 19:00 running: governor up, cap 4; spawning lanes A (amplitude readout), B (hopping Hamiltonian), C (FAIL18 detector + readouts), D (Adversary).
- 2026-09-14 19:05 running: four lanes spawned 19:02 (A, B, C, D). Box: user baseline ~10 GB (chrome 4.5, claude 2.1, Code 1.9); RAM 88 to 91% before any lane job; jobrun will serialise jobs needing more than ~1 GB. Lanes: keep est-ram honest, checkpoint every 10 min, expect waits.
- 2026-09-14 19:10 running: user steer relayed to all four lanes and written as contract addendum 1 (built chain only, no cosmetic variants, D attacks every positive before acceptance).
- next: hourly ledger commits; reassign lanes as they finish; final S28 report in s27/.

## A

## B

## C
- 2026-09-14 19:12 running: contract, brief, S27 report/prereg, router record (S22 L7, S23 L6/L7, S26 L110/L115) read; PREREG_S28_C.md written (no endpoint read). Next: s28_C_fail18.py features + nested logistic + permutation null; s28_C_readout.py; tests; 1-target probe under jobrun.

## D
- 2026-09-14 19:16 running: suite split into TEST jobs (light 8 files 286 pass/3 skip/0 fail, peak 0.91 GB; pipeline + integration queued after two governor kills at 95.7%/98.2% RAM); prereg checks posted S28-L1/L2/L3 (all STANDS WITH CAVEAT); leakage grep on s28_A_amp.py, s28_C_fail18.py, s28_C_readout.py next.
- next: NaN-poison on each lane script as it lands; first S27 reproduction (seed 101, pool row); attack the first positive within the hour.
