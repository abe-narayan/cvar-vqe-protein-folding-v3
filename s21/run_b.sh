#!/bin/sh
# SPRINT 21 / WORKSTREAM B -- one heavy process, BLAS capped at 1, blocks in priority order.
# A and Q first (the coordinator's redirect); the encoding gauge sweep behind them.
cd "$(dirname "$0")/.." || exit 1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
LOG=s21/results/_driver.log
echo "=== driver start $(date) ===" >> "$LOG"
for step in "a:DIST" "q:DIST" "a:LEG" "q:LEG" "a:AMBc" "q:AMBc" "gate" "e:DIST" "en:DIST" "e:AMBc" "en:AMBc"; do
  echo "--- $step $(date) ---" >> "$LOG"
  python -u -m s21.qb3_run "$step" 20 >> "$LOG" 2>&1
  echo "--- $step rc=$? $(date) ---" >> "$LOG"
done
echo "=== driver done $(date) ===" >> "$LOG"
