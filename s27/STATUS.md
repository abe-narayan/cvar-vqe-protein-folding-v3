
# S28 STATUS

## Coordinator
- 2026-09-14 19:00 running: governor up, cap 4; spawning lanes A (amplitude readout), B (hopping Hamiltonian), C (FAIL18 detector + readouts), D (Adversary).
- 2026-09-14 19:05 running: four lanes spawned 19:02 (A, B, C, D). Box: user baseline ~10 GB (chrome 4.5, claude 2.1, Code 1.9); RAM 88 to 91% before any lane job; jobrun will serialise jobs needing more than ~1 GB. Lanes: keep est-ram honest, checkpoint every 10 min, expect waits.
- 2026-09-14 19:10 running: user steer relayed to all four lanes and written as contract addendum 1 (built chain only, no cosmetic variants, D attacks every positive before acceptance).
- next: hourly ledger commits; reassign lanes as they finish; final S28 report in s27/.

## A
- 2026-09-14 19:25 running: prereg + addenda committed (02f1b69b, +1); `s27/s28_A_amp.py`, `s27/s28_A_analyse.py`, `tests/test_s28_A.py` (21 pass); 1-target probe under jobrun OK (3.8 s, 1A13: frame dev 2.1e-14); job `s28A_oracle_126` (ORACLE expressivity ceilings, all 126) running.
- next: ledger S28-L1 (ORACLE ceiling, labelled), then job `s28A_recog_126` (every recognition arm, point cloud), then the built-chain job for the primary list.

## B
- 2026-09-14 19:21 running: `s28B_run` (126 targets, 39 arms each: VQE s0/s1 + eigensolver x {REAL, PERM, RAND} x J {0.1, 0.3, 1, 3} + J = 0; point cloud; resumable JSONL) and `s28B_train` (gradient variance vs J, n = 4..9, 12 targets). Prereg `s27/PREREG_S28_B.md` committed (fee55ba8); module + 14 tests committed (5d81b6d7). J = 0 reproduces S27's DIS VQE arm on 1A13 (2.6124, m = 72). Probe peak RSS 0.341 GB.
- next: post F5 (trainability vs J) and the departure diagnostics to the ledger as soon as `s28B_train` lands; then `--analyse`, chain arms, ledger per arm.

## C
- 2026-09-14 19:12 running: contract, brief, S27 report/prereg, router record (S22 L7, S23 L6/L7, S26 L110/L115) read; PREREG_S28_C.md written (no endpoint read). Next: s28_C_fail18.py features + nested logistic + permutation null; s28_C_readout.py; tests; 1-target probe under jobrun.

## D
- 2026-09-14 19:16 running: suite split into TEST jobs (light 8 files 286 pass/3 skip/0 fail, peak 0.91 GB; pipeline + integration queued after two governor kills at 95.7%/98.2% RAM); prereg checks posted S28-L1/L2/L3 (all STANDS WITH CAVEAT); leakage grep on s28_A_amp.py, s28_C_fail18.py, s28_C_readout.py next.
- next: NaN-poison on each lane script as it lands; first S27 reproduction (seed 101, pool row); attack the first positive within the hour.
