#!/bin/sh
# AGENT C, Sprint 20 -- the remaining heavy work, STRICTLY SEQUENTIAL so this workstream never
# has more than one heavy process on a shared box (BRIEF section 11).
#
# ORDER OF WORK, and why:
#   0. wait for `s20.c_land` (Q2, already running when this driver was re-issued)
#   1. `s20.c_q1relax` -- Q1 RE-SPECIFIED through a shared Relax operator.  Promoted ahead of
#      everything else because the coordinator's LEDGER L6 makes it the Q1 PRIMARY.
#   2. `s20.c_repair` -- AMBER-first vs projection-first (the coordinator's earlier question)
#   3. Q3: merge whatever the other lane's `s19.agentC_pareto` produced and finish the remainder
#
# A SECOND LANE IS RUNNING THE IDENTICAL Q3 EXPERIMENT (`python -m s19.agentC_pareto`, same code,
# same cfg_hash).  I stopped my duplicate rather than theirs and merge their rows in step 3.
#
# Every stage is resumable per target and writes its own `complete` flag, so an interruption
# leaves a named partial and never a silently-truncated result.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$(dirname "$0")/.."

# --- 0. wait for the landscape run already in flight -------------------------------------
while python -c "
import json,sys,os
p='s20/results/c_land.json'
sys.exit(1 if (os.path.exists(p) and json.load(open(p)).get('complete')) else 0)
" 2>/dev/null; do
  sleep 60
done
echo '=== c_land complete ==='
python -u -m s20.c_report land > s20/results/c_land_report.txt 2>&1
python -u -m s20.c_land_null > s20/results/c_land_null_report.txt 2>&1

# --- 1. Q1 re-specified: the shared-operator comparison ----------------------------------
python -u -m s20.c_q1relax > s20/results/c_q1relax.log 2>&1
python -u -m s20.c_q1relax --report > s20/results/c_q1relax_report.txt 2>&1
echo '=== c_q1relax done ==='

# --- 2. the repair path ------------------------------------------------------------------
python -u -m s20.c_repair > s20/results/c_repair.log 2>&1
echo '=== c_repair done ==='

# --- 3. Q3: merge the other lane's rows, compute only the remainder -----------------------
python -u -m s20.c_pareto_finish >> s20/results/c_pareto.log 2>&1
python -u -m s20.c_report pareto > s20/results/c_pareto_report.txt 2>&1
echo '=== ALL DONE ==='
