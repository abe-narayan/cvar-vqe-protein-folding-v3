#!/bin/sh
# AGENT C, Q2.  Gate GC3 first (the declared iteration bound must be INERT for the incumbent
# protocol); only then the n = 126 Pareto run.  Both wait on physical memory: this box is shared
# and `core.amber` refuses to build a Context above 92% physical.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")/.."
python -u -m s19.agentC_pareto --gate > s19/results/agentC_gates_amber.log 2>&1 || exit 1
python -c "import json,sys;o=json.load(open('s19/results/agentC_gates_amber.json'));sys.exit(0 if o['passed'] else 1)" || {
  echo "GATE FAILED -- refusing to run Q2" >> s19/results/agentC_gates_amber.log; exit 1; }
python -u -m s19.agentC_pareto > s19/results/agentC_pareto.log 2>&1
python -u -m s19.agentC_report pareto > s19/results/agentC_pareto_report.txt 2>&1
