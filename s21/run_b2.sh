#!/bin/sh
# SPRINT 21 / WORKSTREAM B -- second-stage driver, re-ordered by scientific priority once
# Blocks A and Q reported on the OpenMM-free objectives.
#
#   1-2  backfill the cross-readout (readout H != training H) on the deployed axis -- cheap
#   3    q:AMBc      F-Q3, the ONE question that needs the badly-conditioned Hamiltonian
#   4-7  the encoding gauge sweep on the OpenMM-free objectives
#   8    the F-E5 harness gate at full n
#   9-11 the AMBER-bound remainder, lowest priority, may be reported PARTIAL
#
# OpenMM contexts are serialised on this box tonight, so every AMBER-bound step is queued behind
# everything that does not need one.  One heavy process, BLAS capped at 1.
cd "$(dirname "$0")/.." || exit 1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
LOG=s21/results/_driver.log
echo "=== driver2 start $(date) ===" >> "$LOG"
for step in "s:DIST" "s:LEG" "q:AMBc" "e:DIST" "en:DIST" "e:LEG" "en:LEG" "gate" "a:DIST" "q:DIST" "s:AMBc" "a:AMBc" "e:AMBc" "en:AMBc"; do
  echo "--- $step $(date) ---" >> "$LOG"
  python -u -m s21.qb3_run "$step" 20 >> "$LOG" 2>&1
  echo "--- $step rc=$? $(date) ---" >> "$LOG"
done
echo "=== driver2 done $(date) ===" >> "$LOG"
