#!/bin/sh
# PHYSICS workstream heavy chain -- ONE process at a time, in priority order.
cd "$(dirname "$0")/.."
set -x
python -m s17.phys_ca --pass 1 >> s17/results/phys_ca_p1.log 2>&1
python -m s17.phys_ca --pass 2 >> s17/results/phys_ca_p2.log 2>&1
python -m s17.phys_ca --pass 4 >> s17/results/phys_ca_p4.log 2>&1
python -m s17.phys_ca --pass 3 >> s17/results/phys_ca_p3.log 2>&1
echo DONE_CHAIN
